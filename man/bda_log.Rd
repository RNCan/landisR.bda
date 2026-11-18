\name{bda_log}
\alias{bda_log}
\title{Combine the Base BDA event summaries}
\usage{bda_log(run, extension = "bda")}
\arguments{
\item{run}{A coupled simulation returned by \code{landisR::landis_sim()}.}
\item{extension}{Name of the Base BDA disturbance in the extension list.}
}
\value{A data frame with one row per epidemic and the Base BDA log fields.}
\description{Reports the regional outbreak status, damaged sites, cohorts
killed and mean severity rounded to two decimal places. Rows are ordered by year
and by agent execution within each year.}

\examples{
\dontrun{
# run is a coupled simulation with an extension named "bda".
outbreaks <- bda_log(run)
head(outbreaks)
}
}
