\name{sampleDataForLtmleMSM}
\alias{sampleDataForLtmleMSM}
\docType{data}
\title{
Sample data, regimes, and summary measures
}
\description{
Sample data for use with ltmleMSM. Data:  n=1000: male age CD4_1 A1 Y1 CD4_2 A2 Y2 CD4_3 A3 Y3 
A1..A3 are treatment nodes, Y1..Y3 are death, CD4_1..CD4_3 are time varying covariates.
We are interested in static regimes where a patient switches at some time.
In summary.measures, switch.time is first time where At is 1 (4 if never switch), time is the horizon.
}
\usage{data(sampleDataForLtmleMSM)}
\format{
  List with three components: data, regimes, summary.measures
}
\details{
	regimes: 200 x 3 x 4 [n x numACnodes x numRegimes]
	summary.measures: 4 x 2 x 3 [numRegimes x numSummaryMeasures x numFinalYnodes]
}
\source{
simulated data
}
\examples{
data(sampleDataForLtmleMSM)
}
\keyword{datasets}
