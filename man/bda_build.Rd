\name{bda_build}
\alias{bda_build}
\title{Build Base BDA parameters from R tables}
\usage{
bda_build(timestep, agents, species, ecoregions = NULL,
  modifiers = NULL, ignored_species = NULL)
}
\arguments{
\item{timestep}{Positive whole number in years, matching the succession step.}
\item{agents}{One row per agent, in execution order. The key column is
  \code{agent}. The other required columns are \code{BDPCalibrator},
  \code{SRDMode}, \code{OutbreakPattern}, \code{TimeSinceLastEpidemic},
  \code{TemporalType}, \code{MinROS}, \code{MaxROS}, \code{Dispersal},
  \code{NeighborFlag}, \code{NeighborRadius}, \code{NeighborWeight},
  \code{IntensityClass2_BDP} and \code{IntensityClass3_BDP}.
  Use \code{MinInterval} and \code{MaxInterval} for \code{CyclicUniform} or
  \code{Mean} and \code{StDev} for \code{CyclicNormal}.
  For spatial dispersal supply \code{DispersalTemplate}, \code{DispersalRate},
  \code{InitialEpicenterNum}, \code{EpidemicThresh},
  \code{OutbreakEpicenterThresh}, \code{SeedEpicenter} and
  \code{SeedEpicenterCoeff}. For neighborhood effects supply
  \code{NeighborShape} and \code{NeighborSpeedUp}. Optional calendar
  columns are \code{StartYear} and \code{EndYear}. In a mixed-agent table
  leave unused conditional settings as \code{NA}. Flags accept logical
  values or the strings \code{"yes"} and \code{"no"}.}
\item{species}{One row per agent and configured species. Columns are
  \code{agent}, \code{species},
  \code{minor_host_age}, \code{minor_host_srd},
  \code{secondary_host_age}, \code{secondary_host_srd},
  \code{primary_host_age}, \code{primary_host_srd},
  \code{resistant_host_age}, \code{resistant_host_vuln},
  \code{tolerant_host_age}, \code{tolerant_host_vuln},
  \code{vulnerable_host_age}, \code{vulnerable_host_vuln} and
  \code{cfs_conifer}. Ages are whole numbers in years. SRD and
  vulnerability coefficients are between zero and one.}
\item{ecoregions}{Optional table with \code{agent}, \code{ecoregion},
  and numeric \code{modifier}. Each agent and ecoregion combination is unique.}
\item{modifiers}{Optional table with \code{agent}, numeric \code{modifier},
  positive integer \code{duration} in years and \code{selectors}. The latter
  is a list-column of character vectors or a column of pipe-delimited
  strings for CSV, for example \code{"FireSeverity3|CUT"}.}
\item{ignored_species}{Optional table with \code{agent} and \code{species}.}
}
\description{Builds a Base BDA parameter set from data frames.
Tables can be created in R or read from CSV. Agent, species and ecoregion
names must match their definitions in the parameter tables and project.}
\value{Base BDA parameters for \code{bda_prepare_coupling()} or
\code{bda_extension()}. The initial calendar consumes random draws in
the row order of \code{agents}.}
\seealso{\code{bda_read()}, \code{bda_prepare_coupling()}}
\examples{
agents <- data.frame(agent = "budworm", BDPCalibrator = 1L, SRDMode = "mean",
  OutbreakPattern = "CyclicUniform", MinInterval = 10, MaxInterval = 20,
  TimeSinceLastEpidemic = 0L, TemporalType = "pulse", MinROS = 1L, MaxROS = 2L,
  Dispersal = FALSE, NeighborFlag = FALSE, NeighborRadius = 100L,
  NeighborWeight = 1, IntensityClass2_BDP = 0.25, IntensityClass3_BDP = 0.50)
species <- data.frame(agent = "budworm", species = "ABIE.BAL",
  minor_host_age = 0L, minor_host_srd = 0.25,
  secondary_host_age = 20L, secondary_host_srd = 0.50,
  primary_host_age = 40L, primary_host_srd = 1,
  resistant_host_age = 11L, resistant_host_vuln = 0.35,
  tolerant_host_age = 20L, tolerant_host_vuln = 0.65,
  vulnerable_host_age = 50L, vulnerable_host_vuln = 1,
  cfs_conifer = TRUE)
parameters <- bda_build(10L, agents, species)
}
