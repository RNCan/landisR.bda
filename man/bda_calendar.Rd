\name{bda_calendar}
\alias{bda_calendar}
\title{Read Base BDA calendar diagnostics}
\usage{bda_calendar(run, extension = "bda")}
\arguments{
\item{run}{A coupled run returned by \code{landisR::landis_sim()}.}
\item{extension}{The name assigned to the Base BDA extension in the
  simulation's list of extensions.}
}
\value{A data frame with one initialization row and one row per evaluated step for
each agent, in execution order. Initial rows appear when the extension is
initialized or first runs. A run without recorded diagnostics
returns an empty data frame.}
\description{Reports the computed interval, elapsed time since the last
pulse, announced next year, pulse flag, regional outbreak status (ROS),
cumulative 32-bit MT19937 word counts immediately before and after the
calendar calculation, and the Boolean cache after the calculation. Set
\code{trace_calendar = TRUE} in \code{bda_extension()} to record these
diagnostics. Multi-agent runs add the agent name and \code{next_year_shared}
as published by the LANDIS-II extension. The \code{next_year} column
reports the next year stored in each agent's own calendar.}
