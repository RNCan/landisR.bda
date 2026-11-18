\name{landisR.bda-package}
\alias{landisR.bda-package}
\alias{landisR.bda}
\title{landisR.bda: Biological Disturbance Agents for LANDIS-R}
\description{
Run Base BDA outbreaks, remove affected cohorts and record outbreak history
on the shared forest before Biomass Succession.
}
\details{
Use \code{\link{bda_read}} to read LANDIS-II input files or
\code{\link{bda_build}} to construct agent and host parameters in R.
Use \code{\link{bda_prepare_coupling}} to initialize the calendar before
other disturbances and register the returned extension in the simulation.
}
\seealso{
\code{\link{bda_extension}}, \code{\link{bda_calendar}},
\code{\link{bda_maps}}, \code{\link{bda_log}}, \code{\link{bda_mortality}},
\code{\link{bda_disturbance_source}}
}
\docType{package}
