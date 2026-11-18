\name{bda_mortality}
\alias{bda_mortality}
\title{Report cohort and biomass mortality from Base BDA}
\usage{bda_mortality(run, extension = "bda", by_species = FALSE)}
\arguments{
\item{run}{A coupled simulation returned by \code{landisR::landis_sim()}.}
\item{extension}{Name of the Base BDA disturbance in the extension list.}
\item{by_species}{If \code{TRUE}, return one row per species killed by each
  epidemic and agent.}
}
\value{A data frame with one row per epidemic, the fields in \code{bda_log()},
and \code{MgBiomassKilled}, the total living biomass killed in megagrams
(metric tonnes). With \code{by_species = TRUE},
columns are \code{Time}, \code{AgentName}, \code{species}, \code{CohortsKilled}
and \code{MgBiomassKilled}. An epidemic without cohort mortality has no species
rows. If the historical species identity is unavailable, the event total is
reported in one row with \code{species = NA}.}
\description{Reports the number of cohorts killed and their biomass immediately
before removal. Rows follow year and agent order. With Biomass Succession,
cohort biomass in grams per square metre is converted to mass using the
project cell area. An epidemic with no deaths has zero killed biomass. A
run without Biomass Succession or a historical event lacking the required
biomass data returns \code{NA} for unavailable biomass.}
\seealso{\code{bda_log()}, \code{bda_prepare_coupling()}}
\examples{
\dontrun{
mortality <- bda_mortality(run, extension = "bda")
write.csv(mortality, "bda_mortality.csv", row.names = FALSE)
mortality_species <- bda_mortality(run, by_species = TRUE)
}
}
