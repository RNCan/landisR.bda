// Synchronous model adapted from LANDIS-II Foundation Base BDA 4.0.1.
// Algorithm references: PlugIn.cs, SiteResources.cs and Epidemic.cs.
#include <Rcpp.h>
#include <landisR/mt19937.hpp>
#include <landisR/progress.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

namespace {
struct Cohorts {
  std::vector<int> cell;
  std::vector<int> species;
  std::vector<double> age;
  std::vector<double> biomass;
  void append(int id, int sp, double a, double b) {
    cell.push_back(id);
    species.push_back(sp + 1);
    age.push_back(a);
    biomass.push_back(b);
  }
  DataFrame table() const {
    return DataFrame::create(_["cell_id"] = cell, _["species"] = species,
                             _["age"] = age, _["biomass"] = biomass);
  }
};

// The v7 succession extension preserves the shared MT state but does not
// carry the boolean cache. Keep that cache in BDA's checkpointed schedule.
struct IntervalRandom {
  landis_rng::Stream& stream;
  std::uint32_t buffer;
  int remaining;

  bool next_boolean() {
    if (remaining == 0) {
      do {
        buffer = static_cast<std::uint32_t>(stream.generator());
        ++stream.draws;
      } while (buffer == std::numeric_limits<std::uint32_t>::max());
      remaining = 31;
    } else {
      buffer >>= 1;
      --remaining;
    }
    return (buffer & 1U) != 0;
  }

  double normal(double mean, double sd) {
    double u, v, radius;
    do {
      u = 2 * stream.next_double() - 1;
      v = 2 * stream.next_double() - 1;
      // Preserve the separate IEEE products used by the C# reference.
      const volatile double u2 = u * u;
      const volatile double v2 = v * v;
      radius = u2 + v2;
    } while (radius >= 1 || radius < 1e-6);
    const double scale = std::sqrt(-2 * std::log(radius) / radius) * sd;
    const volatile double centered = (next_boolean() ? u : v) * scale;
    return centered + mean;
  }
};

std::int64_t next_interval(const List& config, IntervalRandom& rng) {
  const double mean = config["normal_mean"];
  if (!NumericVector::is_na(mean)) {
    const double sd = config.containsElementNamed("normal_sd") ?
      as<double>(config["normal_sd"]) : 0;
    const int timestep = config["timestep"];
    double sample = mean;
    if (sd > 0) {
      rng.normal(mean, sd); // Base BDA discards the first complete sample.
      sample = rng.normal(mean, sd);
    }
    if (!std::isfinite(sample) || sample < std::numeric_limits<int>::min() ||
        sample > std::numeric_limits<int>::max())
      stop("The sampled outbreak interval exceeds the 32-bit calendar range.");
    return std::max<std::int64_t>(0, static_cast<std::int64_t>(sample) - timestep / 2);
  }
  const NumericVector bounds = config["interval"];
  return static_cast<std::int64_t>(bounds[0] + std::trunc(rng.stream.next_double() *
                                                       (bounds[1] - bounds[0])));
}

List calendar_row(const char* phase, int year, std::int64_t interval,
                  std::int64_t elapsed, std::int64_t next_year, bool pulse,
                  int ros, const IntervalRandom& random, std::uint64_t before) {
  return List::create(
    _["phase"] = phase, _["year"] = year, _["interval"] = double(interval),
    _["elapsed"] = double(elapsed), _["next_year"] = double(next_year),
    _["pulse"] = int(pulse), _["ros"] = ros, _["rng_before"] = double(before),
    _["rng_after"] = double(random.stream.draws), _["bit_buffer"] = double(random.buffer),
    _["bit_count"] = random.remaining);
}

struct InitialCalendar {
  std::int64_t interval;
  std::int64_t next_year;
};

InitialCalendar initialize_calendar(const List& config, IntervalRandom& rng,
                                    std::int64_t elapsed) {
  const int start = config["start"];
  const int timestep = config["timestep"];
  const auto interval = next_interval(config, rng) + start;
  return {interval, std::max<std::int64_t>(interval - elapsed, std::max(timestep, start))};
}

int map_percent(double value) {
  if (!std::isfinite(value)) return 0;
  const auto rounded = static_cast<std::int64_t>(std::nearbyint(100 * value));
  // C# stores the map in a signed 16-bit pixel.
  const auto unsigned_code = static_cast<std::uint16_t>(rounded);
  return unsigned_code < 32768 ? unsigned_code : static_cast<int>(unsigned_code) - 65536;
}

enum Zone : unsigned char { no_zone, last_zone, new_zone };

struct SiteDisturbances {
  struct Rule {
    double value;
    int duration;
    CharacterVector labels;
    IntegerVector fire_codes;
  };
  std::vector<Rule> rules;
  std::vector<int> fire_row, harvest_row, bda_row;
  IntegerVector fire_year, fire_severity, harvest_year, bda_year;
  CharacterVector prescriptions, agents;
  String generic_harvest = "Harvest";
  String generic_bda = "BDA";

  SiteDisturbances(const List& config, Nullable<List> history, int n) {
    if (history.isNull() || !config.containsElementNamed("modifiers")) return;
    const List definitions = config["modifiers"];
    if (definitions.size() == 0) return;
    const List sites = history.get();
    if (sites.containsElementNamed("fire")) {
      const List fire = sites["fire"];
      const IntegerVector positions = fire["position"];
      if (positions.size() > 0) {
        fire_year = as<IntegerVector>(fire["year"]);
        fire_severity = as<IntegerVector>(fire["severity"]);
        fire_row.assign(n, -1);
        for (int i = 0; i < positions.size(); ++i) fire_row[positions[i]] = i;
      }
    }
    if (sites.containsElementNamed("harvest")) {
      const List harvest = sites["harvest"];
      const IntegerVector positions = harvest["position"];
      if (positions.size() > 0) {
        harvest_year = as<IntegerVector>(harvest["year"]);
        prescriptions = as<CharacterVector>(harvest["prescription"]);
        harvest_row.assign(n, -1);
        for (int i = 0; i < positions.size(); ++i) harvest_row[positions[i]] = i;
      }
    }
    if (sites.containsElementNamed("bda")) {
      const List bda = sites["bda"];
      const IntegerVector positions = bda["position"];
      if (positions.size() > 0) {
        bda_year = as<IntegerVector>(bda["year"]);
        agents = as<CharacterVector>(bda["agent"]);
        bda_row.assign(n, -1);
        for (int i = 0; i < positions.size(); ++i) bda_row[positions[i]] = i;
      }
    }
    if (fire_row.empty() && harvest_row.empty() && bda_row.empty()) return;
    for (const auto& item : definitions) {
      const List rule = item;
      rules.push_back({as<double>(rule["value"]), as<int>(rule["duration"]),
                       as<CharacterVector>(rule["labels"]), as<IntegerVector>(rule["fire"])});
    }
  }

  double modifier(int position, int year) const {
    const int h = harvest_row.empty() ? -1 : harvest_row[position];
    const int f = fire_row.empty() ? -1 : fire_row[position];
    const int b = bda_row.empty() ? -1 : bda_row[position];
    if (h < 0 && f < 0 && b < 0) return 0;
    const double time_since_harvest = h < 0 ? 0 : double(year) - harvest_year[h];
    const double time_since_fire = f < 0 ? 0 : double(year) - fire_year[f];
    const double time_since_bda = b < 0 ? 0 : double(year) - bda_year[b];
    double total = 0;
    for (const auto& rule : rules) {
      // Preserve input order, including duplicate matches, and the C# order
      // of harvest, fire, and BDA contributions for each row.
      if (h >= 0 && time_since_harvest <= rule.duration) {
        const double contribution = rule.value * std::max(0.0, time_since_harvest) / rule.duration;
        const String prescription = prescriptions[h];
        for (const String name : rule.labels)
          if (name == prescription || name == generic_harvest) total += contribution;
      }
      if (f >= 0 && time_since_fire <= rule.duration) {
        const double contribution = rule.value * std::max(0.0, time_since_fire) / rule.duration;
        for (int code : rule.fire_codes)
          if (code == -1 || code == fire_severity[f]) total += contribution;
      }
      if (b >= 0 && time_since_bda > 0 && time_since_bda <= rule.duration) {
        const double contribution = rule.value * time_since_bda / rule.duration;
        const String agent = agents[b];
        for (const String name : rule.labels)
          if (name == agent || name == generic_bda) total += contribution;
      }
    }
    return total;
  }
};

struct Grid {
  int rows, cols;
  double cell_length;
  IntegerVector index, ids;

  int site(int row, int col) const {
    if (row < 0 || row >= rows || col < 0 || col >= cols) return -1;
    return index[row + col * rows];
  }
  int row(int i) const { return (ids[i] - 1) % rows; }
  int col(int i) const { return (ids[i] - 1) / rows; }
  int sample(landis_rng::Stream& rng) const {
    const int sampled_row = rng.next(rows);
    const int sampled_col = rng.next(cols);
    return site(sampled_row, sampled_col);
  }
};

List disperse(const List& settings, const List& schedule, const Grid& grid,
              const IntegerVector& order, const IntegerMatrix& deltas,
              const std::vector<double>& vulnerability, double distance,
              landis_rng::Stream& rng, std::vector<Zone>& zone, bool trace) {
  const IntegerVector previous = schedule.containsElementNamed("previous_zone") ?
    as<IntegerVector>(schedule["previous_zone"]) : IntegerVector(0);
  for (int i : previous) zone[i] = last_zone;
  const bool first = previous.size() == 0;
  const int old_count = schedule.containsElementNamed("epicenter_count") ?
    as<int>(schedule["epicenter_count"]) : as<int>(settings["initial_count"]);
  const auto before = rng.draws;
  int inside = 0;
  int outside = 0;
  bool can_seed = false;
  const double threshold = settings["threshold"];
  if (!first) {
    const IntegerVector damaged = schedule["previous_damage"];
    const double inside_threshold = settings["inside_threshold"];
    for (int i : damaged)
      if (vulnerability[i] >= inside_threshold) ++inside;
    for (int i : order) {
      if (zone[i] != last_zone && vulnerability[i] >= threshold) {
        ++outside;
        can_seed |= vulnerability[i] > threshold;
      }
    }
  }
  // In 4.0.1 shuffle empties its input and the caller discards its return.
  // The inside candidates consume draws but none are selected.
  for (int remaining = inside; remaining > 0; --remaining) rng.next(remaining);
  const auto after_inside = rng.draws;
  std::vector<int> centers;
  if (!first && as<bool>(settings["seed"]) && outside > 0) {
    const double coefficient = settings["seed_coefficient"];
    int remaining = static_cast<int>(outside *
      (coefficient <= 0 ? 1 : std::exp(-coefficient * old_count)));
    if (remaining > 0 && (!can_seed || coefficient <= 0))
      stop("Base BDA external epicenter selection cannot finish with these thresholds or SeedEpicenterCoeff.");
    while (remaining > 0) {
      const int i = grid.sample(rng);
      if (i < 0 || zone[i] != no_zone || !(vulnerability[i] > threshold)) continue;
      centers.push_back(i);
      remaining = static_cast<int>(remaining * std::exp(-coefficient * centers.size()));
    }
  }
  const int outside_selected = centers.size();
  const auto after_outside = rng.draws;
  if (centers.empty()) {
    for (int k = 0; k < old_count;) {
      const int i = grid.sample(rng);
      if (i < 0) continue;
      centers.push_back(i);
      ++k;
    }
  }
  const auto after_selection = rng.draws;
  const bool fixed_radius = as<std::string>(settings["template"]) == "MaxRadius";
  std::vector<int> queue;
  for (int center : centers) {
    const int center_row = grid.row(center);
    const int center_col = grid.col(center);
    if (fixed_radius) {
      for (int k = 0; k < deltas.nrow(); ++k) {
        const int i = grid.site(center_row + deltas(k, 0), center_col + deltas(k, 1));
        if (i >= 0) zone[i] = new_zone;
      }
      continue;
    }
    queue.clear();
    queue.push_back(center);
    zone[center] = new_zone;
    for (std::size_t k = 0; k < queue.size(); ++k) {
      const int row = grid.row(queue[k]);
      const int col = grid.col(queue[k]);
      for (int j = 0; j < deltas.nrow(); ++j) {
        const int r = row + deltas(j, 0);
        const int c = col + deltas(j, 1);
        const int i = grid.site(r, c);
        if (i < 0 || zone[i] == new_zone || !(vulnerability[i] > threshold)) continue;
        const double dr = r - center_row;
        const double dc = c - center_col;
        if (std::sqrt(dr * dr + dc * dc) * grid.cell_length > distance) continue;
        queue.push_back(i);
        zone[i] = new_zone;
      }
    }
  }
  List result = List::create(_["epicenter_count"] = static_cast<int>(centers.size()));
  if (trace) {
    IntegerVector center_cells(centers.size());
    CharacterVector sources(centers.size(), "fallback");
    for (int k = 0; k < center_cells.size(); ++k) {
      center_cells[k] = grid.ids[centers[k]];
      if (k < outside_selected) sources[k] = "outside";
    }
    result["trace"] = List::create(
      _["first"] = int(first), _["old_count"] = old_count, _["inside_candidates"] = inside,
      _["outside_candidates"] = outside, _["outside_selected"] = outside_selected,
      _["selected"] = static_cast<int>(centers.size()), _["rng_before"] = double(before),
      _["rng_after_inside"] = double(after_inside), _["rng_after_outside"] = double(after_outside),
      _["rng_after_selection"] = double(after_selection),
      _["epicenters"] = DataFrame::create(_["cell_id"] = center_cells, _["source"] = sources));
  }
  return result;
}
}

// [[Rcpp::export(rng = false)]]
List bda_prime_cpp(const List& config, const List& rng_state, bool trace_calendar,
                   double bit_buffer = 0, int bit_count = 0, int year = 0) {
  landis_rng::Stream rng(0);
  rng.draws = landis_rng::Stream::offset(as<double>(rng_state["offset"]));
  rng.restore(as<std::string>(rng_state["state"]));
  const auto before = rng.draws;
  IntervalRandom random{rng, static_cast<std::uint32_t>(bit_buffer), bit_count};
  const auto elapsed = static_cast<std::int64_t>(as<double>(config["initial_elapsed"]));
  const auto initial = initialize_calendar(config, random, elapsed);
  RObject trace = R_NilValue;
  if (trace_calendar)
    trace = calendar_row("initial", year, initial.interval, elapsed, initial.next_year,
                         false, 0, random, before);
  return List::create(
    _["schedule"] = List::create(
      _["initialized"] = true, _["interval"] = double(initial.interval),
      _["elapsed"] = double(elapsed), _["next_year"] = double(initial.next_year),
      _["bit_buffer"] = double(random.buffer), _["bit_count"] = random.remaining),
    _["rng"] = List::create(_["state"] = rng.encode(), _["offset"] = double(rng.draws)),
    _["calendar"] = trace);
}

// [[Rcpp::export(rng = false)]]
List bda_step_cpp(const List& forest, const List& config, const List& schedule,
                  const List& topology, const List& rng_state, int year, bool trace_calendar = false,
                  bool trace_dispersal = false, Nullable<List> disturbances = R_NilValue) {
  landis_progress::Reporter progress;
  progress.phase(0, "Outbreak calendar");
  landis_rng::Stream rng(0);
  rng.draws = landis_rng::Stream::offset(as<double>(rng_state["offset"]));
  rng.restore(as<std::string>(rng_state["state"]));
  const auto rng_before = rng.draws;
  const int timestep = config["timestep"];
  auto elapsed = static_cast<std::int64_t>(as<double>(schedule["elapsed"]));
  auto interval = static_cast<std::int64_t>(as<double>(schedule["interval"]));
  auto next_year = static_cast<std::int64_t>(as<double>(schedule["next_year"]));

  IntervalRandom calendar_rng{
    rng,
    schedule.containsElementNamed("bit_buffer") ?
      static_cast<std::uint32_t>(as<double>(schedule["bit_buffer"])) : 0U,
    schedule.containsElementNamed("bit_count") ? as<int>(schedule["bit_count"]) : 0
  };
  List calendar;
  if (!as<bool>(schedule["initialized"])) {
    const auto initial = initialize_calendar(config, calendar_rng, elapsed);
    interval = initial.interval;
    next_year = initial.next_year;
    if (trace_calendar)
      calendar.push_back(calendar_row("initial", year - timestep, interval, elapsed,
                                      next_year, false, 0, calendar_rng, rng_before));
  }
  const auto rng_before_check = rng.draws;
  elapsed += timestep;
  if (elapsed > 10000)
    stop("TimeSinceLastEpidemic is outside the LANDIS-II range [0, 10000].");
  int ros = config["min_ros"];
  const bool pulse = interval <= elapsed && year <= as<int>(config["end"]);
  if (pulse) {
    interval = next_interval(config, calendar_rng);
    next_year = year + interval;
    elapsed = 0;
    const int max_ros = config["max_ros"];
    if (as<bool>(config["variable"]))
      ros += 1 + static_cast<int>(rng.next_double() * (max_ros - ros));
    else
      ros = max_ros;
  }
  const auto rng_after_calendar = rng.draws;
  if (trace_calendar)
    calendar.push_back(calendar_row("step", year, interval, elapsed, next_year, pulse,
                                    ros, calendar_rng, rng_before_check));
  const bool spatial = config.containsElementNamed("dispersal") &&
    !Rf_isNull(config["dispersal"]);
  List next_schedule = List::create(
    _["initialized"] = true, _["interval"] = double(interval), _["elapsed"] = double(elapsed),
    _["next_year"] = double(next_year), _["bit_buffer"] = double(calendar_rng.buffer),
    _["bit_count"] = calendar_rng.remaining);
  if (spatial && schedule.containsElementNamed("previous_zone")) {
    next_schedule["previous_zone"] = schedule["previous_zone"];
    next_schedule["previous_damage"] = schedule["previous_damage"];
    next_schedule["epicenter_count"] = schedule["epicenter_count"];
  }
  List output = List::create(
    _["schedule"] = next_schedule,
    _["calendar"] = calendar,
    _["pulse"] = pulse,
    _["event"] = R_NilValue,
    _["rng"] = List::create(_["state"] = rng_after_calendar == rng_before ?
                            as<std::string>(rng_state["state"]) : rng.encode(),
                            _["offset"] = double(rng_after_calendar)));
  if (ros <= 0) return output;

  const List states = forest["states"];
  const IntegerVector cell_state = forest["cell_state"];
  const IntegerVector ids = forest["cell_id"];
  const IntegerVector offset = states["offset"];
  const IntegerVector species = states["species_id"];
  const NumericVector ages = states["age"];
  const NumericVector biomass = states["biomass"];
  const NumericMatrix parameters = topology["parameters"];
  const LogicalVector ignored = topology["ignored"];
  const NumericVector eco = topology["eco"];
  const IntegerVector order = topology["order"];
  const LogicalVector trace_index = topology["trace"];
  const int n = ids.size();
  const int nrow = forest["nrow"];
  const int ncol = forest["ncol"];
  const int base = offset[0];
  const Grid grid{nrow, ncol, as<double>(forest["cell_length"]),
                  as<IntegerVector>(topology["index"]), ids};
  const SiteDisturbances site_disturbances(config, disturbances, n);

  // Compute the resource once for each interned CSR state that is used.
  std::vector<double> by_state(offset.size() - 1, 0);
  std::vector<bool> computed(offset.size() - 1, false);
  std::vector<double> oldest(parameters.nrow(), 0);
  std::vector<double> srd(n), modified(n), nrd(n, 0), vulnerability(n);
  const bool mean_srd = config["mean_srd"];
  progress.start(0.08, 0.30, "Host resources (SRD)", n);
  for (int i = 0; i < n; ++i) {
    progress.advance(i);
    const int state = cell_state[i] - 1;
    if (!computed[state]) {
      std::fill(oldest.begin(), oldest.end(), 0);
      for (int j = offset[state] - base; j < offset[state+1] - base; ++j)
        oldest[species[j]] = std::max(oldest[species[j]], ages[j]);
      double sum = 0;
      double max = 0;
      int present = 0;
      for (int sp = 0; sp < parameters.nrow(); ++sp) {
        if (oldest[sp] <= 0 || ignored[sp]) continue;
        ++present;
        double host = 0;
        for (int grade = 0; grade < 3; ++grade)
          if (oldest[sp] >= parameters(sp, grade)) host = parameters(sp, grade + 6);
        sum += host;
        max = std::max(max, host);
      }
      by_state[state] = mean_srd ? (present ? sum / present : R_NaN) : max;
      computed[state] = true;
    }
    srd[i] = by_state[state];
    modified[i] = srd[i] > 0 ?
      std::clamp(srd[i] + site_disturbances.modifier(i, year) + eco[i], 0.0, 1.0) : 0;
  }

  const double scale = (double(ros) / 3) * as<double>(config["calibrator"]);
  const double weight = config["neighbor_weight"];
  std::vector<Zone> zone(spatial ? n : 0, no_zone);
  List dispersion;
  if (spatial) {
    progress.phase(0.30, "Outbreak dispersal");
    for (int i = 0; i < n; ++i) vulnerability[i] = std::max(0.0, scale * modified[i]);
    const List settings = config["dispersal"];
    const double distance = as<double>(settings["rate"]) * timestep;
    dispersion = disperse(settings, schedule, grid, order, topology["dispersal_deltas"],
                          vulnerability, distance, rng, zone, trace_dispersal);
  }
  if (as<bool>(config["neighbor"])) {
    const IntegerMatrix deltas = topology["deltas"];
    const NumericVector resource_weights = topology.containsElementNamed("resource_weights") ?
      as<NumericVector>(topology["resource_weights"]) : NumericVector(deltas.nrow(), 1.0);
    const int stride = config.containsElementNamed("neighbor_stride") ?
      as<int>(config["neighbor_stride"]) : 1;
    struct Neighbor { int site, offset; };
    std::vector<Neighbor> neighborhood;
    neighborhood.reserve(deltas.nrow());
    progress.start(0.47, 0.77, "Neighborhood (NRD)", order.size());
    std::size_t visited = 0;
    for (int i : order) {
      progress.advance(visited++);
      if (!(srd[i] > 0)) continue;
      neighborhood.clear();
      if (!spatial || zone[i] == new_zone) {
        const int row = grid.row(i);
        const int col = grid.col(i);
        for (int j = 0; j < deltas.nrow(); ++j) {
          const int position = grid.site(row + deltas(j, 0), col + deltas(j, 1));
          if (position >= 0) neighborhood.push_back({position, j});
        }
      }
      const bool keep_all = neighborhood.size() <= static_cast<std::size_t>(stride);
      double total = 0, weight_sum = 0;
      int shuffled_position = 0;
      // C# shuffles every active neighbor before applying the stride.
      while (!neighborhood.empty()) {
        const int position = rng.next(static_cast<int>(neighborhood.size()));
        if (keep_all || shuffled_position % stride == 0) {
          const Neighbor neighbor = neighborhood[position];
          const double neighbor_weight = resource_weights[neighbor.offset];
          const volatile double weighted_resource = modified[neighbor.site] * neighbor_weight;
          total += weighted_resource;
          weight_sum += neighbor_weight;
        }
        neighborhood.erase(neighborhood.begin() + position);
        ++shuffled_position;
      }
      nrd[i] = weight_sum > 0 ? total / weight_sum : 0;
      // LANDIS-II rounds the weighted NRD before adding the site resource.
      const volatile double weighted_nrd = nrd[i] * weight;
      vulnerability[i] = std::max(0.0, scale * ((modified[i] + weighted_nrd) / (1 + weight)));
    }
  } else if (!spatial) {
    for (int i = 0; i < n; ++i) vulnerability[i] = std::max(0.0, scale * modified[i]);
  }
  const auto rng_after_neighbors = rng.draws;
  const NumericVector thresholds = config["thresholds"];
  Cohorts survivors, killed;
  std::vector<int> damaged, severity, dead_count, cfs_count;
  std::vector<double> dead_biomass;
  std::vector<int> trace_ids, map_severity, map_srd, map_nrd, map_bdp;
  std::vector<int> trace_zone, next_zone, next_damage;
  std::vector<double> trace_srd, trace_modified, trace_nrd, trace_bdp, trace_draw;
  std::vector<bool> death;
  double severity_total = 0;
  progress.start(0.77, 0.95, "Severity & mortality", order.size());
  std::size_t visited = 0;
  for (int i : order) {
    progress.advance(visited++);
    const double u = rng.next_double();
    const double bdp = vulnerability[i];
    int final_severity = 0;
    if (spatial && zone[i] == new_zone) next_zone.push_back(i);
    if ((!spatial || zone[i] == new_zone) && bdp > u) {
      const int state = cell_state[i] - 1;
      const int first = offset[state] - base;
      const int last = offset[state+1] - base;
      int deaths = 0;
      int cfs = 0;
      double removed = 0;
      death.assign(last - first, false);
      for (int j = first; j < last; ++j) {
        const int sp = species[j];
        for (int grade = 0; grade < 3; ++grade) {
          if (ages[j] >= parameters(sp, grade+3) && u <= bdp * parameters(sp, grade+9)) {
            death[j-first] = true;
            break;
          }
        }
        if (death[j-first]) {
          ++deaths;
          cfs += parameters(sp, 12) != 0;
          removed += biomass[j];
          killed.append(ids[i], sp, ages[j], biomass[j]);
        }
      }
      if (deaths) {
        if (spatial) next_damage.push_back(i);
        final_severity = bdp >= thresholds[1] ? 3 : bdp >= thresholds[0] ? 2 : 1;
        severity_total += final_severity;
        damaged.push_back(ids[i]);
        severity.push_back(final_severity);
        dead_count.push_back(deaths);
        cfs_count.push_back(cfs);
        dead_biomass.push_back(removed);
        for (int j = first; j < last; ++j)
          if (!death[j-first]) survivors.append(ids[i], species[j], ages[j], biomass[j]);
      }
    }
    if (trace_index.size() && trace_index[i]) {
      trace_ids.push_back(ids[i]);
      if (spatial) trace_zone.push_back(zone[i]);
      map_severity.push_back(1 + final_severity);
      map_srd.push_back(map_percent(srd[i]));
      map_nrd.push_back(map_percent(nrd[i]));
      map_bdp.push_back(map_percent(bdp));
      trace_srd.push_back(srd[i]);
      trace_modified.push_back(modified[i]);
      trace_nrd.push_back(nrd[i]);
      trace_bdp.push_back(bdp);
      trace_draw.push_back(u);
    }
  }
  progress.phase(0.96, "Writing outbreak diagnostics");
  output["event"] = List::create(
    _["ros"] = ros, _["mean_severity"] = damaged.empty() ? 0 : severity_total / damaged.size(),
    _["sites"] = DataFrame::create(_["cell_id"] = damaged, _["severity"] = severity,
                                    _["cohorts_killed"] = dead_count, _["cfs_killed"] = cfs_count,
                                    _["biomass_killed"] = dead_biomass),
    _["cohorts_killed"] = killed.table(),
    _["offsets"] = NumericVector::create(_["before"] = double(rng_before),
                                        _["after_calendar"] = double(rng_after_calendar),
                                        _["after_neighbors"] = double(rng_after_neighbors),
                                        _["after_sites"] = double(rng.draws)),
    _["trace"] = DataFrame::create(_["cell_id"] = trace_ids, _["severity_code"] = map_severity,
                                  _["srd_code"] = map_srd, _["nrd_code"] = map_nrd,
                                  _["bdp_code"] = map_bdp, _["srd"] = trace_srd,
                                  _["modified_srd"] = trace_modified, _["nrd"] = trace_nrd,
                                  _["bdp"] = trace_bdp, _["uniform"] = trace_draw));
  output["survivors"] = survivors.table();
  if (spatial) {
    next_schedule["previous_zone"] = wrap(next_zone);
    next_schedule["previous_damage"] = wrap(next_damage);
    next_schedule["epicenter_count"] = dispersion["epicenter_count"];
    output["schedule"] = next_schedule;
    List event = output["event"];
    DataFrame trace_sites = event["trace"];
    trace_sites.push_back(trace_zone, "zone");
    event["trace"] = trace_sites;
    if (trace_dispersal) event["dispersion"] = dispersion["trace"];
    output["event"] = event;
  }
  output["rng"] = List::create(_["state"] = rng.encode(), _["offset"] = double(rng.draws));
  return output;
}
