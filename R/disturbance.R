bda_publish_disturbance <- function(simulation, type, sites) {
  type <- match.arg(type, c("fire", "harvest"))
  forest <- landisR::landis_csr(simulation)
  quantity <- if (type == "fire") "severity" else "cohorts_killed"
  columns <- c("cell_id", "year", quantity, if (type == "harvest") "prescription")
  if (!is.data.frame(sites) || !all(columns %in% names(sites)))
    stop("sites requires columns ", paste(columns, collapse = ", "), ".", call. = FALSE)
  position <- match(sites$cell_id, forest$cell_id)
  if (anyNA(position) || anyDuplicated(position))
    stop("sites must contain distinct active cell IDs.", call. = FALSE)
  is_integer <- function(x) {
    is.numeric(x) && all(is.finite(x) & abs(x) <= .Machine$integer.max & trunc(x) == x)
  }
  values <- sites[[quantity]]
  if (!is_integer(sites$year) || !is_integer(values) || any(values < 0))
    stop("Disturbance years and nonnegative severity or cohort counts must be integers.",
         call. = FALSE)
  record <- list(position = position[values > 0] - 1L,
                 year = as.integer(sites$year[values > 0]))
  if (type == "fire") {
    if (any(values > 255L)) stop("Fire severity must be in [0, 255].", call. = FALSE)
    record$severity <- as.integer(values[values > 0])
  } else {
    if (!is.character(sites$prescription) || anyNA(sites$prescription))
      stop("Harvest prescription names must be strings.", call. = FALSE)
    record$prescription <- enc2utf8(trimws(sites$prescription[values > 0]))
    record$cohorts_killed <- as.integer(values[values > 0])
  }
  simulation$bda_disturbances[[type]] <- record
  simulation
}

bda_disturbance_source <- function(extension, type) {
  type <- match.arg(type, c("fire", "harvest"))
  if (!inherits(extension, "landisR_extension") || extension$type != "disturbance")
    stop("extension must be a landisR disturbance extension.", call. = FALSE)
  force(extension)
  step <- function(simulation, state, year) {
    result <- extension$step(simulation, state, year)
    events <- result$events
    if (is.null(events)) return(result)
    if (type == "fire") {
      sites <- events$burned
      if (!is.data.frame(sites) || !all(c("cell_id", "severity") %in% names(sites)))
        stop("The fire extension must return its burned cell table.", call. = FALSE)
    } else {
      if (!is.data.frame(events) || !all(c("cells", "Prescription") %in% names(events)))
        stop("The harvest extension must return its cell selections and prescriptions.",
             call. = FALSE)
      cells <- as.integer(unlist(events$cells, use.names = FALSE))
      prescriptions <- rep(as.character(events$Prescription), lengths(events$cells))
      if (length(cells)) {
        before <- landisR::landis_csr(simulation)
        after <- landisR::landis_csr(result$simulation)
        position <- match(cells, before$cell_id)
        if (anyNA(position) || anyDuplicated(position))
          stop("Harvest selections must refer to distinct active cells.", call. = FALSE)
        cohort_count <- function(forest) {
          index <- forest$cell_state[position]
          forest$states$offset[index + 1L] - forest$states$offset[index]
        }
        killed <- cohort_count(before) - cohort_count(after)
        # Biomass Harvest schedules establishment for succession, so a
        # partial removal keeps the cohort and a complete removal drops it.
        keep <- killed > 0L
        cells <- cells[keep]
        prescriptions <- prescriptions[keep]
        killed <- killed[keep]
      } else killed <- integer()
      sites <- data.frame(cell_id = cells, cohorts_killed = killed,
                          prescription = prescriptions)
    }
    sites$year <- rep.int(as.integer(year), nrow(sites))
    result$simulation <- bda_publish_disturbance(result$simulation, type, sites)
    result
  }
  landisR::landis_extension(step, extension$initialize, type = "disturbance")
}
