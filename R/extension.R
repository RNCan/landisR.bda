bda_dispersal_deltas <- function(template, distance, cell_length, nrow, ncol) {
  if (template == "MaxRadius") {
    reach <- floor(distance / cell_length)
    dy <- min(reach, nrow - 1L)
    dx <- min(reach, ncol - 1L)
    grid <- expand.grid(col = seq.int(-dx, dx), row = seq.int(-dy, dy))
    within <- sqrt((grid$row * cell_length)^2 + (grid$col * cell_length)^2) <= distance
    deltas <- as.matrix(grid[within, c("row", "col")])
  } else {
    # The templates are cumulative in the order used by GetNeighbors().
    deltas <- rbind(c(0, 1), c(1, 0), c(0, -1), c(-1, 0),
                    c(-1, 1), c(1, 1), c(1, -1), c(-1, -1),
                    c(-2, 0), c(0, 2), c(2, 0), c(0, -2),
                    c(-2, -2), c(-2, -1), c(-1, -2), c(-2, 2), c(-2, 1), c(-1, 2),
                    c(2, 2), c(1, 2), c(2, 1), c(2, -2), c(2, -1), c(1, -2))
    deltas <- deltas[seq_len(as.integer(sub("N$", "", template))), , drop = FALSE]
  }
  storage.mode(deltas) <- "integer"
  deltas
}

bda_check_initial_simulation <- function(simulation, parameters) {
  if (simulation$year != 0L || is.null(simulation$landis_rng))
    stop("Initialize at year zero with recruitment_rng = 'landis_v7'.", call. = FALSE)
  if (simulation$project$canonical$timestep != parameters$timestep)
    stop("Base BDA and Biomass Succession must use the same timestep.", call. = FALSE)
}

bda_grid <- function(forest, trace_cells) {
  ids <- forest$cell_id
  row <- (ids - 1L) %% forest$nrow + 1L
  col <- (ids - 1L) %/% forest$nrow + 1L
  index <- rep.int(-1L, forest$nrow * forest$ncol)
  index[ids] <- seq_along(ids) - 1L
  if (identical(trace_cells, "all")) trace_cells <- ids
  list(order = order(row, col, method = "radix") - 1L, index = index,
       trace = if (length(trace_cells)) ids %in% trace_cells else logical())
}

bda_agent_state <- function(simulation, parameters, forest) {
  nrow_grid <- forest$nrow
  ncol_grid <- forest$ncol
  labels <- simulation$project$canonical$landscape$cells$ecoregion[forest$cell_id]
  configured_species <- rownames(parameters$species)
  unknown <- setdiff(c(configured_species, parameters$ignored), forest$species)
  if (length(unknown))
    stop("Unknown species in Base BDA: ", paste(unknown, collapse = ", "), call. = FALSE)
  eco_labels <- simulation$project$canonical$ecoregions$ecoregion
  unknown_eco <- setdiff(names(parameters$eco), eco_labels)
  if (length(unknown_eco))
    stop("Unknown Base BDA ecoregion: ", paste(unknown_eco, collapse = ", "), call. = FALSE)
  eco <- unname(parameters$eco[match(labels, names(parameters$eco))])
  eco[is.na(eco)] <- 0

  species <- matrix(rep(c(rep(999, 6L), rep(0, 7L)), each = length(forest$species)),
                    nrow = length(forest$species))
  species[match(configured_species, forest$species), ] <- parameters$species
  resource_weights <- numeric()
  if (parameters$neighbor) {
    # LANDIS-II stores cell length and Gaussian weights as single precision.
    cell_length <- bda_single(forest$cell_length)
    reach <- floor(parameters$radius / cell_length)
    dy <- min(reach, nrow_grid - 1L)
    dx <- min(reach, ncol_grid - 1L)
    grid <- expand.grid(col = seq.int(-dx, dx), row = seq.int(-dy, dy))
    distance <- sqrt((grid$col * cell_length)^2 + (grid$row * cell_length)^2)
    within <- distance > 0 & distance <= parameters$radius
    deltas <- as.matrix(grid[within, c("row", "col")])
    storage.mode(deltas) <- "integer"
    distance <- distance[within]
    shape <- if (is.null(parameters$neighbor_shape)) "uniform" else parameters$neighbor_shape
    resource_weights <- switch(shape,
      uniform = rep(1, length(distance)),
      linear = 1 - distance / parameters$radius,
      gaussian = bda_single(exp(-(distance^2) / (parameters$radius %/% 2)^2)))
  } else deltas <- matrix(integer(), ncol = 2L)
  dispersal <- parameters$dispersal
  dispersal_deltas <- NULL
  if (!is.null(dispersal)) {
    distance <- as.double(dispersal$rate) * parameters$timestep
    if (distance > .Machine$integer.max)
      stop("DispersalRate times Timestep exceeds the 32-bit distance range.", call. = FALSE)
    dispersal_deltas <- bda_dispersal_deltas(dispersal$template, distance, forest$cell_length,
                                             nrow_grid, ncol_grid)
  }
  list(schedule = list(initialized = FALSE, elapsed = parameters$initial_elapsed,
                       interval = 0, next_year = 0),
       topology = list(parameters = species, ignored = forest$species %in% parameters$ignored,
                       eco = eco, deltas = deltas, resource_weights = resource_weights,
                       dispersal_deltas = dispersal_deltas))
}

bda_initialize <- function(simulation, parameters, trace_cells) {
  bda_check_initial_simulation(simulation, parameters)
  forest <- landisR::landis_csr(simulation)
  state <- bda_agent_state(simulation, parameters, forest)
  topology <- c(bda_grid(forest, trace_cells), state$topology)
  state$topology <- topology[c("parameters", "ignored", "eco", "order", "index", "deltas",
                               "resource_weights", "dispersal_deltas", "trace")]
  state
}

bda_single <- function(values) {
  readBin(writeBin(values, raw(), size = 4L), "double", n = length(values), size = 4L)
}

bda_prepare_coupling <- function(simulation, parameters, trace_cells = integer(),
                                 trace_calendar = FALSE, trace_dispersal = FALSE) {
  extension <- bda_extension(parameters, trace_cells, trace_calendar, trace_dispersal)
  bda_check_initial_simulation(simulation, parameters)
  initial <- if (is.null(parameters$agents))
    bda_prime_cpp(parameters, simulation$landis_rng, trace_calendar)
  else
    bda_prime_agents(parameters$agents, simulation$landis_rng, trace_calendar)
  simulation$landis_rng <- initial$rng
  extension$initialize <- bda_prepared_initializer(extension$initialize, initial)
  list(simulation = simulation, extension = extension)
}

bda_prepared_initializer <- function(initialize, initial) {
  force(initialize)
  force(initial)
  function(simulation) {
    if (simulation$year != 0L || !identical(simulation$landis_rng, initial$rng))
      stop("Use the simulation returned by bda_prepare_coupling().", call. = FALSE)
    state <- initialize(simulation)
    if (!is.null(state$agents)) return(bda_adopt_prime(state, initial, 0L))
    state$schedule <- initial$schedule
    if (!is.null(initial$calendar)) state$calendar <- list("0" = list(initial$calendar))
    state
  }
}

bda_extension <- function(parameters, trace_cells = integer(), trace_calendar = FALSE,
                          trace_dispersal = FALSE) {
  if (!inherits(parameters, "landisR_bda_parameters"))
    stop("Use bda_read() or bda_build() to create Base BDA parameters.", call. = FALSE)
  if (!identical(trace_cells, "all") &&
      (!is.numeric(trace_cells) || anyNA(trace_cells) ||
       any(!is.finite(trace_cells) | trace_cells < 1 | trace_cells != trunc(trace_cells))))
    stop("trace_cells must be positive cell IDs or 'all'.", call. = FALSE)
  if (!identical(trace_calendar, TRUE) && !identical(trace_calendar, FALSE))
    stop("trace_calendar must be TRUE or FALSE.", call. = FALSE)
  if (!identical(trace_dispersal, TRUE) && !identical(trace_dispersal, FALSE))
    stop("trace_dispersal must be TRUE or FALSE.", call. = FALSE)
  force(parameters)
  force(trace_cells)
  force(trace_calendar)
  force(trace_dispersal)
  if (!is.null(parameters$agents))
    return(bda_multi_extension(parameters, trace_cells, trace_calendar, trace_dispersal))

  history <- bda_initial_history(list(parameters), parameters$agent)
  initialize <- function(simulation) {
    state <- bda_initialize(simulation, parameters, trace_cells)
    state$history <- history
    state
  }
  step <- function(simulation, state, year) {
    result <- bda_run_agent(simulation, parameters, state$schedule, state$topology, year,
                            trace_calendar, trace_dispersal, state$history)
    state$schedule <- result$schedule
    if (trace_calendar) state$calendar[[as.character(year)]] <- result$calendar
    state$history <- bda_update_history(state$history, result$event, state$topology$index)
    list(simulation = result$simulation, state = state, events = result$event)
  }
  structure(landisR::landis_extension(step, initialize, type = "disturbance"),
            progress_label = "Base BDA")
}

bda_run_agent <- function(simulation, parameters, schedule, topology, year,
                          trace_calendar, trace_dispersal, history,
                          report = getOption("landisR.progress.report")) {
  progress <- if (!is.null(report)) function(fraction, label)
    report(fraction, paste0(parameters$agent, ": ", label))
  if (!is.null(progress)) progress(0.01, "Reading hosts")
  forest <- landisR::landis_csr(simulation)
  disturbances <- simulation$bda_disturbances
  if (!is.null(history) && length(history$position)) disturbances$bda <- history
  result <- landisR::landis_progress_native(progress, 0.05, 0.84,
    bda_step_cpp(forest, parameters, schedule, topology, simulation$landis_rng,
                 year, trace_calendar, trace_dispersal, disturbances))
  event <- result$event
  simulation$landis_rng <- result$rng
  if (!is.null(event)) {
    event$year <- as.integer(year)
    event$agent <- parameters$agent
    event$cohorts_killed$species <- forest$species[event$cohorts_killed$species]
    event$dimensions <- c(forest$nrow, forest$ncol)
    event$active_cells <- length(forest$cell_id)
    if (nrow(event$sites)) {
      if (!is.null(progress)) progress(0.86, "Applying cohort mortality")
      survivors <- result$survivors
      survivors$species <- forest$species[survivors$species]
      simulation <- landisR::landis_apply_cohort_batch(
        simulation, survivors, event$sites$cell_id,
        event_metadata = list(list(extension = "Base BDA", year = year,
                                   agent = parameters$agent, cells = event$sites$cell_id)),
        preserve_cohort_order = TRUE)
      if (!is.null(progress)) progress(0.93, "Scheduling regeneration")
      simulation <- landisR::landis_schedule_regeneration(simulation, event$cohorts_killed)
    }
  }
  if (!is.null(progress)) progress(0.98, "Updating outbreak history")
  list(simulation = simulation, schedule = result$schedule, event = event,
       calendar = result$calendar, pulse = result$pulse)
}

bda_maps <- function(event) {
  if (is.null(event$trace) || nrow(event$trace) != event$active_cells)
    stop("Create the extension with trace_cells = 'all' to reconstruct maps.", call. = FALSE)
  fields <- c(severity = "severity_code", srd = "srd_code",
              nrd = "nrd_code", bdp = "bdp_code")
  result <- lapply(fields, function(field) {
    values <- matrix(0L, nrow = event$dimensions[1L], ncol = event$dimensions[2L])
    values[event$trace$cell_id] <- event$trace[[field]]
    values
  })
  result
}

bda_log_events <- function(run, extension) {
  if (!inherits(run, "landisR_coupled_run") || !extension %in% names(run$events))
    stop("Expected a landis_sim() run containing the named Base BDA extension.",
         call. = FALSE)
  events <- unname(run$events[[extension]])
  if (!is.null(run$extension_states[[extension]]$agents))
    events <- do.call(c, events)
  events
}

bda_format_log <- function(events) {
  data.frame(
    Time = vapply(events, function(e) e$year, integer(1L)),
    ROS = vapply(events, function(e) e$ros, integer(1L)),
    AgentName = vapply(events, function(e) e$agent, character(1L)),
    CohortsKilled = vapply(events, function(e) nrow(e$cohorts_killed), integer(1L)),
    DamagedSites = vapply(events, function(e) nrow(e$sites), integer(1L)),
    MeanSeverity = vapply(events, function(e) floor(100 * e$mean_severity + 0.5) / 100,
                          numeric(1L))
  )
}

bda_log <- function(run, extension = "bda") {
  bda_format_log(bda_log_events(run, extension))
}

bda_killed_density <- function(event) {
  if (!nrow(event$cohorts_killed)) return(0)
  density <- event$sites$biomass_killed
  if (is.null(density)) density <- event$cohorts_killed$biomass
  if (is.null(density)) NA_real_ else sum(density)
}

bda_mortality_by_species <- function(events, area_mg) {
  empty <- data.frame(Time = integer(), AgentName = character(), species = character(),
                      CohortsKilled = integer(), MgBiomassKilled = numeric())
  rows <- lapply(events, function(event) {
    killed <- event$cohorts_killed
    n <- nrow(killed)
    if (!n) return(NULL)
    if (is.null(killed$species) || anyNA(killed$species)) {
      return(data.frame(Time = event$year, AgentName = event$agent, species = NA_character_,
                        CohortsKilled = n, MgBiomassKilled = bda_killed_density(event) * area_mg))
    }
    density <- if (is.null(killed$biomass)) rep(NA_real_, n) else killed$biomass
    groups <- rowsum(cbind(cohorts = rep.int(1L, n), density = density),
                     killed$species)
    data.frame(Time = event$year, AgentName = event$agent, species = rownames(groups),
               CohortsKilled = as.integer(groups[, "cohorts"]),
               MgBiomassKilled = unname(groups[, "density"]) * area_mg)
  })
  result <- do.call(rbind, c(list(empty), rows))
  rownames(result) <- NULL
  result
}

bda_mortality <- function(run, extension = "bda", by_species = FALSE) {
  if (!isTRUE(by_species) && !identical(by_species, FALSE))
    stop("by_species must be TRUE or FALSE.", call. = FALSE)
  events <- bda_log_events(run, extension)
  simulation <- run$simulation
  # Biomass Succession stores cohort density in g/m2. Return Mg per epidemic.
  area_mg <- if (identical(simulation$growth_engine$id, "biomass52"))
    simulation$project$canonical$cell_length^2 / 1e6 else NA_real_
  if (by_species) return(bda_mortality_by_species(events, area_mg))
  log <- bda_format_log(events)
  log$MgBiomassKilled <- vapply(events, bda_killed_density, numeric(1L)) * area_mg
  log
}

bda_calendar <- function(run, extension = "bda") {
  if (!inherits(run, "landisR_coupled_run") || !extension %in% names(run$extension_states))
    stop("Expected a landis_sim() run containing the named Base BDA extension.",
         call. = FALSE)
  entries <- run$extension_states[[extension]]$calendar
  if (!length(entries)) return(data.frame())
  rows <- do.call(c, unname(entries))
  do.call(rbind, lapply(rows, as.data.frame))
}
