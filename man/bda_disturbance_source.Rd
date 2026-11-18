\name{bda_disturbance_source}
\alias{bda_disturbance_source}
\alias{bda_publish_disturbance}
\title{Publish disturbance site variables to Base BDA}
\usage{
bda_disturbance_source(extension, type)
bda_publish_disturbance(simulation, type, sites)
}
\arguments{
\item{extension}{A \code{landisR} disturbance extension. The fire adapter
requires the Base Fire burned-cell table. The harvest adapter requires the
Biomass Harvest cell selections and prescription names.}
\item{type}{Either \code{"fire"} or \code{"harvest"}.}
\item{simulation}{A \code{landisR} simulation.}
\item{sites}{A data frame with unique active \code{cell_id} values and
integer event \code{year} values. Include \code{severity} for fire or
\code{cohorts_killed} and character \code{prescription} for harvest.}
}
\description{
Pass the adapted extension to \code{landisR::landis_sim()} in its original
position relative to Base BDA. For another producer of site variables, use
\code{bda_publish_disturbance()} directly from its disturbance step.
}
\details{
Publishes only positive severity or complete-removal counts. A publication
replaces the previous report of the same type, including when the report is
empty. A skipped producer step leaves the existing report available to BDA.
The simulation checkpoint stores the current sparse reports.

The optional Harvest adapter bridges complete cohort removals to the
\code{Harvest.CohortsKilled} field that Base BDA 4.0.1 reads. With stock
Biomass Harvest 4.4, the native field is named \code{Harvest.CohortsDamaged}.
Use the plain Harvest extension to reproduce the stock field availability.
No additional RNG draws are consumed by these adapters.
}
\value{A wrapped disturbance extension or an updated simulation.}

\examples{
grid <- data.frame(cell_id = 1L, row = 1L, col = 1L, active = TRUE)
project <- landisR::landis_build_landscape(grid, duration = 10L,
  cell_length = 100)
sim <- landisR::landis_initialize(project)
sim <- bda_publish_disturbance(sim, type = "fire",
  sites = data.frame(cell_id = 1L, year = 0L, severity = 3L))

# A disturbance can also publish the burned-cell table it returns.
producer <- landisR::landis_extension(function(simulation, state, year) {
  burned <- data.frame(cell_id = 1L, severity = 2L)
  list(simulation = simulation, state = state, events = list(burned = burned))
})
fire_source <- bda_disturbance_source(producer, type = "fire")
run <- landisR::landis_sim(sim, extensions = list(fire = fire_source),
  progress = FALSE)
}
