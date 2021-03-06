% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_stics.R
\name{run_stics}
\alias{run_stics}
\title{Run STICS}
\usage{
run_stics(dirpath = getwd(), exe_name = "stics")
}
\arguments{
\item{dirpath}{USM directory path}

\item{exe_name}{(optional) The executable name.}
}
\value{
The function prints the STICS output to the console and returns
\code{TRUE} if STICS ran successfully, or an error if any problem occured
}
\description{
Helper function to start a simulation using
STICS executable
}
\examples{
\dontrun{
library(sticRs)
run_stics()
}
}
