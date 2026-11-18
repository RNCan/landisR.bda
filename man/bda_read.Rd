\name{bda_read}
\alias{bda_read}
\title{Read LANDIS-II Base BDA input}
\usage{bda_read(path, scenario_directory = dirname(path))}
\arguments{
\item{path}{Path to the \code{Base BDA} input file.}
\item{scenario_directory}{Working directory of the LANDIS-II scenario. Relative
  agent paths are resolved from this location.}
}
\value{Base BDA parameters for \code{bda_extension()}.}
\description{Reads one or more agents for Base BDA 4.0.1. The port covers
synchronous and spatial outbreaks, uniform and normal intervals, start and end years,
pulse or variable pulse, mean or maximum site resource dominance, ecoregion
modifiers, and uniform, linear, or Gaussian resource neighborhoods.
Use \code{NeighborSpeedUp} with \code{none}, \code{2x}, \code{3x}, or
\code{4x} to select how neighbors contribute after shuffling. Spatial dispersal
supports the \code{MaxRadius}, \code{4N}, \code{8N}, \code{12N}, and
\code{24N} templates. The \code{DisturbanceModifiers} table accepts
an effect, a positive integer duration, and one or more selectors per row.
Fire selectors include \code{Fire} and \code{FireSeverity3}. Harvest
selectors include \code{Harvest} and case-sensitive prescription names.
Use \code{bda_disturbance_source()} to publish those events. BDA selectors
match \code{BDA} or a case-sensitive agent name. They refer to the latest
BDA mortality event on each site.
}
\details{Put the first agent file on the \code{BDAInputFiles} line. Put each
additional file on its own line immediately after it. Agent names must be
distinct. The file order governs initialization and execution within each step.}
