bda_tokens <- function(path) {
  lines <- trimws(sub("<<.*$", "", readLines(path, warn = FALSE)))
  lines <- lines[nzchar(lines) & !startsWith(lines, ">>")]
  lapply(lines, function(line) scan(text = line, what = character(), quiet = TRUE,
                                   quote = '"', comment.char = ""))
}

bda_value <- function(rows, key, default = NULL) {
  result <- Filter(function(row) identical(row[[1L]], key), rows)
  if (!length(result) && !is.null(default)) return(default)
  if (length(result) != 1L || length(result[[1L]]) != 2L)
    stop("Expected exactly one value for ", key, call. = FALSE)
  result[[1L]][[2L]]
}

bda_number <- function(x, key, lower = -Inf, upper = Inf, integer = FALSE) {
  value <- suppressWarnings(as.numeric(x))
  if (any(!is.finite(value) | value < lower | value > upper) ||
      (integer && any(value != trunc(value) | abs(value) > .Machine$integer.max)))
    stop("Invalid numeric value for ", key, call. = FALSE)
  if (integer) as.integer(value) else value
}

bda_flag <- function(x, key) {
  if (!x %in% c("yes", "no"))
    stop(key, " must be yes or no.", call. = FALSE)
  x == "yes"
}

bda_read_agent <- function(path) bda_agent_from_rows(bda_tokens(path))

bda_agent_from_rows <- function(rows) {
  keys <- vapply(rows, function(row) row[[1L]], character(1))
  locate <- function(key) match(key, keys, nomatch = length(keys) + 1L)
  section <- function(left, right) {
    if (right <= left + 1L) list() else rows[seq.int(left + 1L, right - 1L)]
  }
  if (locate("IntensityClass3_BDP") > length(rows))
    stop("Missing IntensityClass3_BDP.", call. = FALSE)
  table_index <- locate("BDASpeciesParameters")
  if (table_index > length(rows)) stop("Missing BDASpeciesParameters.", call. = FALSE)
  ignored_index <- locate("IgnoredSpecies")
  modifier_index <- locate("DisturbanceModifiers")
  eco_rows <- section(locate("IntensityClass3_BDP"), min(table_index, modifier_index))
  if (any(lengths(eco_rows) != 2L))
    stop("Ecoregion modifiers must contain a name and a value.", call. = FALSE)
  eco_names <- vapply(eco_rows, `[[`, character(1), 1L)
  if (anyDuplicated(eco_names)) stop("Duplicate ecoregion modifier.", call. = FALSE)
  eco <- structure(
    bda_number(vapply(eco_rows, `[[`, character(1), 2L), "ecoregion modifiers", -1, 1),
    names = eco_names
  )
  modifiers <- lapply(section(modifier_index, table_index), function(row) {
    if (length(row) < 3L)
      stop("DisturbanceModifiers rows need a modifier, duration and at least one selector.",
           call. = FALSE)
    names <- row[-c(1L, 2L)]
    labels <- enc2utf8(trimws(names))
    specific <- startsWith(names, "FireSeverity")
    suffix <- substr(names, nchar(names), nchar(names))
    fire <- integer(length(names))
    fire[!specific & labels == "Fire"] <- -1L
    matched <- specific & suffix %in% as.character(1:9)
    fire[matched] <- as.integer(suffix[matched])
    list(value = bda_number(row[1L], "DisturbanceModifier"),
         duration = bda_number(row[2L], "Duration", lower = 1, integer = TRUE),
         labels = labels, fire = fire)
  })
  species_rows <- section(table_index, ignored_index)
  if (any(lengths(species_rows) != 14L))
    stop("Each species row needs a name, 12 numeric values and a CFS flag.", call. = FALSE)
  species <- vapply(species_rows, `[[`, character(1), 1L)
  if (anyDuplicated(species)) stop("Duplicate BDASpeciesParameters species.", call. = FALSE)
  values <- t(vapply(species_rows, function(row) {
    c(bda_number(row[c(2L, 4L, 6L, 8L, 10L, 12L)], "species ages", 0, 999, TRUE),
      bda_number(row[c(3L, 5L, 7L, 9L, 11L, 13L)], "species probabilities", 0, 1),
      bda_flag(row[[14L]], "CFS flag"))
  }, numeric(13L)))
  rownames(values) <- species
  ignored <- section(ignored_index, length(rows) + 1L)
  if (any(lengths(ignored) != 1L))
    stop("Each IgnoredSpecies row must contain one species.", call. = FALSE)
  ignored <- unlist(ignored, use.names = FALSE)

  scalar <- rows[seq_len(locate("IntensityClass3_BDP"))]
  val <- function(key, default = NULL) bda_value(scalar, key, default)
  number <- function(key, default = NULL, lower = 0, upper = Inf, integer = FALSE) {
    bda_number(val(key, default), key, lower, upper, integer)
  }
  if (val("LandisData") != "BDA Agent") stop("Expected LandisData BDA Agent.", call. = FALSE)
  dispersal <- NULL
  if (bda_flag(val("Dispersal"), "Dispersal")) {
    template <- val("DispersalTemplate")
    if (!template %in% c("MaxRadius", "4N", "8N", "12N", "24N"))
      stop("Invalid DispersalTemplate.", call. = FALSE)
    dispersal <- list(template = template, rate = number("DispersalRate", lower = 1, integer = TRUE),
                      initial_count = number("InitialEpicenterNum", integer = TRUE),
                      threshold = number("EpidemicThresh", upper = 1),
                      inside_threshold = number("OutbreakEpicenterThresh", upper = 1),
                      seed = bda_flag(val("SeedEpicenter"), "SeedEpicenter"),
                      seed_coefficient = number("SeedEpicenterCoeff", lower = -Inf))
  }
  neighbor <- bda_flag(val("NeighborFlag"), "NeighborFlag")
  neighbor_shape <- if (neighbor) val("NeighborShape") else "uniform"
  neighbor_stride <- if (neighbor)
    match(val("NeighborSpeedUp"), c("none", "2x", "3x", "4x")) else 1L
  if (!neighbor_shape %in% c("uniform", "linear", "gaussian") || is.na(neighbor_stride))
    stop("Invalid NeighborShape or NeighborSpeedUp.", call. = FALSE)
  srd_mode <- val("SRDMode")
  temporal <- val("TemporalType")
  if (!srd_mode %in% c("mean", "max") || !temporal %in% c("pulse", "variablepulse"))
    stop("Invalid SRDMode or TemporalType.", call. = FALSE)
  pattern <- val("OutbreakPattern")
  if (pattern == "CyclicUniform") {
    interval <- round(c(number("MinInterval"), number("MaxInterval")))
    if (interval[1L] > interval[2L] || any(interval > .Machine$integer.max))
      stop("Invalid CyclicUniform bounds.", call. = FALSE)
    normal_mean <- NA_real_
    normal_sd <- 0
  } else if (pattern == "CyclicNormal") {
    normal_sd <- number("StDev")
    normal_mean <- number("Mean", lower = -.Machine$integer.max, upper = .Machine$integer.max)
    interval <- c(0, 0)
  } else {
    stop("Supported outbreak patterns are CyclicUniform and CyclicNormal.",
         call. = FALSE)
  }
  ros <- c(number("MinROS", integer = TRUE, upper = 3),
           number("MaxROS", integer = TRUE, upper = 3))
  if (ros[1L] > ros[2L]) stop("MinROS exceeds MaxROS.", call. = FALSE)
  list(agent = val("BDAAgentName"), mean_srd = srd_mode == "mean",
       calibrator = number("BDPCalibrator", integer = TRUE),
       start = number("StartYear", "0", integer = TRUE),
       end = number("EndYear", as.character(.Machine$integer.max), integer = TRUE),
       initial_elapsed = number("TimeSinceLastEpidemic", upper = 10000, integer = TRUE),
       interval = interval, normal_mean = normal_mean, normal_sd = normal_sd,
       variable = temporal == "variablepulse", dispersal = dispersal,
       min_ros = ros[1L], max_ros = ros[2L], neighbor = neighbor,
       radius = number("NeighborRadius", integer = TRUE),
       neighbor_shape = neighbor_shape, neighbor_stride = neighbor_stride,
       neighbor_weight = number("NeighborWeight"),
       thresholds = c(number("IntensityClass2_BDP"), number("IntensityClass3_BDP")),
       species = values, ignored = ignored, eco = eco, modifiers = modifiers)
}

bda_read <- function(path, scenario_directory = dirname(path)) {
  rows <- bda_tokens(path)
  if (bda_value(rows, "LandisData") != "Base BDA")
    stop("Expected LandisData Base BDA.", call. = FALSE)
  keys <- vapply(rows, function(row) row[[1L]], character(1))
  i <- match("BDAInputFiles", keys)
  if (is.na(i) || length(rows[[i]]) != 2L || any(lengths(rows[-seq_len(i)]) != 1L))
    stop("BDAInputFiles must name the first agent file, followed by one file per line.",
         call. = FALSE)
  refs <- c(rows[[i]][[2L]], unlist(rows[-seq_len(i)], use.names = FALSE))
  refs <- chartr("\\", "/", refs)
  relative <- !grepl("^(/|[A-Za-z]:/)", refs)
  refs[relative] <- file.path(scenario_directory, refs[relative])
  timestep <- bda_number(bda_value(rows[seq_len(i)], "Timestep"), "Timestep", 1, integer = TRUE)
  bda_parameter_set(lapply(refs, bda_read_agent), timestep)
}

bda_parameter_set <- function(agents, timestep) {
  agents <- lapply(agents, function(agent) {
    agent$timestep <- timestep
    class(agent) <- "landisR_bda_parameters"
    agent
  })
  if (length(agents) == 1L) return(agents[[1L]])
  names(agents) <- vapply(agents, `[[`, character(1L), "agent")
  if (anyDuplicated(names(agents))) stop("Base BDA agent names must be distinct.", call. = FALSE)
  structure(list(timestep = timestep, agents = agents),
            class = c("landisR_bda_multi_parameters", "landisR_bda_parameters"))
}
