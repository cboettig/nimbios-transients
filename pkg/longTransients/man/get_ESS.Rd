% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.R
\name{get_ESS}
\alias{get_ESS}
\title{Compute summaries of effective sample size for derivative of potential for parametric model}
\usage{
get_ESS(samples, x = seq(0.1, 1.9, length.out = 100), na.rm = TRUE)
}
\arguments{
\item{samples}{}

\item{x}{}

\item{na.rm}{default `TRUE`}
}
\value{
vector of min, mean, and median ESS across derivative of potential
}
\description{
Compute summaries of effective sample size for derivative of potential for parametric model
}
