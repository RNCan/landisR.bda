bda_initial_history <- function(agents, agent_names) {
  selectors <- c("BDA", agent_names)
  for (agent in agents) {
    if (any(vapply(agent$modifiers, function(rule) any(rule$labels %in% selectors),
                   logical(1L))))
      return(list(position = integer(), year = integer(), agent = character()))
  }
  NULL
}

bda_update_history <- function(history, event, index) {
  if (is.null(history) || is.null(event) || !nrow(event$sites)) return(history)
  positions <- index[event$sites$cell_id]
  keep <- !history$position %in% positions
  list(position = c(history$position[keep], positions),
       year = c(history$year[keep], rep.int(event$year, length(positions))),
       agent = c(history$agent[keep], rep.int(event$agent, length(positions))))
}

bda_prime_agents <- function(agents, rng, trace_calendar, year = 0L) {
  cache <- list(bit_buffer = 0, bit_count = 0L)
  schedules <- calendar <- vector("list", length(agents))
  for (i in seq_along(agents)) {
    initial <- bda_prime_cpp(agents[[i]], rng, trace_calendar, cache$bit_buffer,
                             cache$bit_count, year)
    rng <- initial$rng
    cache <- initial$schedule[names(cache)]
    schedules[[i]] <- initial$schedule[!names(initial$schedule) %in% names(cache)]
    if (trace_calendar) {
      row <- initial$calendar
      row$agent <- names(agents)[i]
      row$next_year_shared <- row$next_year
      calendar[[i]] <- row
    }
  }
  list(schedules = schedules, cache = cache, rng = rng, calendar = calendar,
       next_year_shared = schedules[[length(agents)]]$next_year)
}

bda_adopt_prime <- function(state, initial, year) {
  for (i in seq_along(state$agents)) state$agents[[i]]$schedule <- initial$schedules[[i]]
  state$cache <- initial$cache
  state$next_year_shared <- initial$next_year_shared
  if (any(lengths(initial$calendar))) state$calendar[[as.character(year)]] <- initial$calendar
  state
}

bda_multi_extension <- function(parameters, trace_cells, trace_calendar, trace_dispersal) {
  agents <- parameters$agents
  history <- bda_initial_history(agents, names(agents))
  force(trace_cells)
  force(trace_calendar)
  force(trace_dispersal)

  initialize <- function(simulation) {
    bda_check_initial_simulation(simulation, parameters)
    forest <- landisR::landis_csr(simulation)
    list(agents = lapply(agents, function(agent) bda_agent_state(simulation, agent, forest)),
         grid = bda_grid(forest, trace_cells), history = history)
  }

  step <- function(simulation, state, year) {
    report <- getOption("landisR.progress.report")
    if (is.null(state$cache)) {
      initial_year <- year - parameters$timestep
      initial <- bda_prime_agents(agents, simulation$landis_rng, trace_calendar, initial_year)
      simulation$landis_rng <- initial$rng
      state <- bda_adopt_prime(state, initial, initial_year)
    }
    events <- calendar <- structure(vector("list", length(agents)), names = names(agents))
    sd_fields <- names(state$cache)
    for (i in seq_along(agents)) {
      agent_report <- if (!is.null(report)) function(fraction, label)
        report((i - 1 + fraction) / length(agents), label)
      agent_state <- state$agents[[i]]
      result <- bda_run_agent(simulation, agents[[i]], c(agent_state$schedule, state$cache),
                              c(state$grid, agent_state$topology), year,
                              trace_calendar, trace_dispersal, state$history, agent_report)
      simulation <- result$simulation
      state$cache <- result$schedule[sd_fields]
      state$agents[[i]]$schedule <- result$schedule[!names(result$schedule) %in% sd_fields]
      state$history <- bda_update_history(state$history, result$event, state$grid$index)
      if (result$pulse) state$next_year_shared <- result$schedule$next_year
      if (trace_calendar) {
        row <- result$calendar[[1L]]
        row$agent <- names(agents)[i]
        row$next_year_shared <- state$next_year_shared
        calendar[[i]] <- row
      }
      events[i] <- list(result$event)
    }
    if (trace_calendar) state$calendar[[as.character(year)]] <- unname(calendar)
    events <- Filter(Negate(is.null), events)
    list(simulation = simulation, state = state,
         events = if (length(events)) events else NULL)
  }

  structure(landisR::landis_extension(step, initialize, type = "disturbance"),
            progress_label = "Base BDA")
}
