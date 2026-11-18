\name{bda_extension}
\alias{bda_extension}
\alias{bda_maps}
\title{Run Base BDA before Biomass Succession}
\usage{
bda_extension(parameters, trace_cells = integer(), trace_calendar = FALSE,
  trace_dispersal = FALSE)
bda_maps(event)
}
\arguments{
\item{parameters}{Parameters returned by \code{bda_read()} or \code{bda_build()}.}
\item{trace_cells}{LANDIS-R cell IDs for which to keep diagnostics at each
  outbreak, or \code{"all"} to enable map reconstruction. The default keeps
  only damaged sites and removed cohorts.}
\item{trace_calendar}{Record the interval, regional outbreak status (ROS)
  and RNG diagnostics at initialization and at each simulation step.
  Retrieve the table with \code{bda_calendar()}.}
\item{trace_dispersal}{With spatial dispersal, record ordered epicenter
  cell IDs and sources, candidate counts, and MT19937 word counts in the
  event's \code{dispersion} list. A separate \code{trace_cells} setting
  records zone membership (0 outside, 1 previous only, 2 reached) for
  the selected cells in the event's \code{trace} data frame.}
\item{event}{One Base BDA outbreak event stored by
  \code{landisR::landis_sim()}. Set \code{trace_cells = "all"} before
  running the extension to obtain the data needed for map reconstruction.}
}
\value{\code{bda_extension()} returns an external disturbance extension.
\code{bda_maps()} returns a list of four integer matrices for severity, SRD,
NRD, and BDP, with rows and columns matching the original grid.}
\description{Applies whole-cohort removal to the shared forest and advances
the same random stream as Biomass Succession. Set
\code{recruitment_rng = "landis_v7"} when initializing the project. Use equal
timesteps for Base BDA and Biomass Succession.}
\details{Use \code{bda_prepare_coupling()} to run the initial calendar draws
before the first disturbance step. A multi-agent run stores events for each
year by agent name, for example \code{run$events$bda[["10"]][["budworm"]]}.
Single-agent runs keep an event directly under each year.
Each agent reads cohort mortality applied by earlier agents in that step.
The coupled progress bar identifies Base BDA while it runs. Control
display with the progress argument of \code{landisR::landis_sim()}.
Call \code{bda_log()} for the combined event table.}

\examples{
\dontrun{
# Use an initialized biomass simulation in "landis_v7" mode.
extension <- bda_extension(parameters, trace_cells = "all")
run <- landisR::landis_sim(sim, extensions = list(bda = extension))

# For a single agent, select a year with an outbreak using bda_log().
outbreaks <- bda_log(run)
if (nrow(outbreaks)) {
  year <- as.character(outbreaks$Time[1L])
  maps <- bda_maps(run$events$bda[[year]])
  maps$severity
}
}
}
