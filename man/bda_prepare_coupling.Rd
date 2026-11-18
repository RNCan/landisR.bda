\name{bda_prepare_coupling}
\alias{bda_prepare_coupling}
\title{Initialize the Base BDA Calendar Before Disturbances}
\usage{
bda_prepare_coupling(simulation, parameters, trace_cells = integer(),
  trace_calendar = FALSE, trace_dispersal = FALSE)
}
\arguments{
  \item{simulation}{A year-zero simulation initialized by
    \code{landisR::landis_initialize()} with \code{recruitment_rng = "landis_v7"}.}
  \item{parameters}{Parameters returned by \code{bda_read()} or \code{bda_build()}.}
  \item{trace_cells, trace_calendar, trace_dispersal}{Optional diagnostics as
    described for \code{bda_extension()}.}
}
\value{A list containing the updated \code{simulation} and the prepared
\code{extension}.}
\description{Draws the first outbreak interval for each agent in file order at year
zero, before an earlier-running disturbance can consume random words.
The calendar and shared Boolean cache then travel with the extension state.}
\details{Use the returned simulation and extension together in
\code{landisR::landis_sim()}. Place the extension at the desired position in the
disturbance list. Resume with the saved simulation and extension state
through the regular checkpoint API.}
