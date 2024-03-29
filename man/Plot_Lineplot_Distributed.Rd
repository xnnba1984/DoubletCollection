% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualization.R
\name{Plot_Lineplot_Distributed}
\alias{Plot_Lineplot_Distributed}
\title{Lineplots for distributed computing}
\usage{
Plot_Lineplot_Distributed(result, data, measurement)
}
\arguments{
\item{result}{A data frame with approriate format to visualize.}

\item{data}{Name of the dataset to visualize.}

\item{measurement}{"AUPRC" or "AUROC".}
}
\value{
No return values.
}
\description{
Plot lineplots of detection accuracy across different batch numbers and methods.
}
\examples{
Plot_Lineplot_Distributed(result = table.batch, data = 'pbmc-ch', measurement = 'AUPRC')

}
