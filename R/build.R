bda_build_table <- function(table, required, label, agents = NULL) {
  if (!is.data.frame(table) || anyDuplicated(names(table)) ||
      !all(required %in% names(table)))
    stop(label, " must be a data.frame with columns: ", paste(required, collapse = ", "),
         call. = FALSE)
  scalar <- setdiff(required, "selectors")
  if (!all(vapply(table[scalar], function(column) is.atomic(column) && is.null(dim(column)),
                  logical(1L))))
    stop(label, " columns must contain scalar values.", call. = FALSE)
  if (!is.null(agents) && any(!table$agent %in% agents))
    stop(label, " contains an unknown agent.", call. = FALSE)
}

bda_build_text <- function(x) {
  if (is.logical(x)) return(ifelse(x, "yes", "no"))
  if (is.numeric(x)) return(sprintf("%.17g", x))
  as.character(x)
}

bda_build <- function(timestep, agents, species, ecoregions = NULL,
                      modifiers = NULL, ignored_species = NULL) {
  if (length(timestep) != 1L)
    stop("Timestep must contain exactly one value.", call. = FALSE)
  timestep <- bda_number(bda_build_text(timestep), "Timestep", 1, integer = TRUE)
  required <- c("agent", "BDPCalibrator", "SRDMode", "OutbreakPattern",
                "TimeSinceLastEpidemic", "TemporalType", "MinROS", "MaxROS",
                "Dispersal", "NeighborFlag", "NeighborRadius", "NeighborWeight",
                "IntensityClass2_BDP", "IntensityClass3_BDP")
  optional <- c("StartYear", "EndYear", "MinInterval", "MaxInterval", "Mean", "StDev",
                "DispersalTemplate", "DispersalRate", "InitialEpicenterNum",
                "EpidemicThresh", "OutbreakEpicenterThresh", "SeedEpicenter",
                "SeedEpicenterCoeff", "NeighborShape", "NeighborSpeedUp")
  bda_build_table(agents, required, "agents")
  if (!nrow(agents) || !is.character(agents$agent) || anyNA(agents$agent) ||
      any(!nzchar(trimws(agents$agent))) || anyDuplicated(agents$agent))
    stop("agents must contain one row per distinct nonempty agent name.", call. = FALSE)
  unknown <- setdiff(names(agents), c(required, optional))
  if (length(unknown))
    stop("Unknown agents columns: ", paste(unknown, collapse = ", "), call. = FALSE)
  if (!all(vapply(agents, function(column) is.atomic(column) && is.null(dim(column)),
                  logical(1L))))
    stop("agents columns must contain scalar values.", call. = FALSE)

  species_columns <- c("species", "minor_host_age", "minor_host_srd",
    "secondary_host_age", "secondary_host_srd", "primary_host_age", "primary_host_srd",
    "resistant_host_age", "resistant_host_vuln", "tolerant_host_age", "tolerant_host_vuln",
    "vulnerable_host_age", "vulnerable_host_vuln", "cfs_conifer")
  bda_build_table(species, c("agent", species_columns), "species", agents$agent)
  if (!is.null(ecoregions))
    bda_build_table(ecoregions, c("agent", "ecoregion", "modifier"), "ecoregions", agents$agent)
  if (!is.null(modifiers))
    bda_build_table(modifiers, c("agent", "modifier", "duration", "selectors"), "modifiers", agents$agent)
  if (!is.null(ignored_species))
    bda_build_table(ignored_species, c("agent", "species"), "ignored_species", agents$agent)

  params <- lapply(seq_len(nrow(agents)), function(i) {
    name <- agents$agent[i]
    settings <- as.list(agents[i, setdiff(names(agents), "agent"), drop = FALSE])
    settings <- settings[!vapply(settings, is.na, logical(1L))]
    missing <- setdiff(required[-1L], names(settings))
    if (length(missing))
      stop("Missing settings for ", name, ": ", paste(missing, collapse = ", "), call. = FALSE)
    scalar_columns <- setdiff(names(settings), "IntensityClass3_BDP")
    scalar <- Map(function(key, value) c(key, bda_build_text(value)),
                  scalar_columns, settings[scalar_columns])
    rows <- c(list(c("LandisData", "BDA Agent"), c("BDAAgentName", name)),
              unname(scalar), list(c("IntensityClass3_BDP",
                                    bda_build_text(settings$IntensityClass3_BDP))))

    eco <- ecoregions[which(ecoregions$agent == name), , drop = FALSE]
    eco_rows <- Map(function(key, value) c(as.character(key), bda_build_text(value)),
                    eco$ecoregion, eco$modifier)
    rows <- c(rows, unname(eco_rows))

    effects <- modifiers[which(modifiers$agent == name), , drop = FALSE]
    if (NROW(effects)) {
      selectors <- effects$selectors
      if (is.factor(selectors)) selectors <- as.character(selectors)
      if (is.character(selectors)) selectors <- strsplit(selectors, "|", fixed = TRUE)
      if (!is.list(selectors) || length(selectors) != nrow(effects) ||
          !all(vapply(selectors, function(x) is.character(x) && length(x) > 0L &&
                      !anyNA(x) && all(nzchar(trimws(x))), logical(1L))))
        stop("modifiers$selectors must contain nonempty selector lists.", call. = FALSE)
      effect_rows <- Map(function(value, duration, names)
        c(bda_build_text(value), bda_build_text(duration), trimws(names)),
        effects$modifier, effects$duration, selectors)
      rows <- c(rows, list("DisturbanceModifiers"), unname(effect_rows))
    }

    cohorts <- species[species$agent == name, species_columns, drop = FALSE]
    species_rows <- lapply(seq_len(nrow(cohorts)), function(j)
      vapply(cohorts[j, , drop = FALSE], bda_build_text, character(1L), USE.NAMES = FALSE))
    ignored <- as.character(ignored_species$species[ignored_species$agent == name])
    if (anyNA(cohorts$species) || any(!nzchar(trimws(cohorts$species))) ||
        anyNA(eco$ecoregion) || any(!nzchar(trimws(eco$ecoregion))) ||
        anyNA(ignored) || any(!nzchar(trimws(ignored))))
      stop("Species and ecoregion names must be nonempty.", call. = FALSE)
    rows <- c(rows, list("BDASpeciesParameters"), species_rows,
              list("IgnoredSpecies"), as.list(ignored))
    bda_agent_from_rows(rows)
  })
  bda_parameter_set(params, timestep)
}
