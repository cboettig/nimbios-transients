% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{simulate_model_nM}
\alias{simulate_model_nM}
\title{Simulate population trajectories}
\usage{
simulate_model_nM(constants = NULL, inits = NULL, INCLUDE_ME = T, seed = 270)
}
\arguments{
\item{constants}{list with constants for nimbleModel()}

\item{inits}{list with initial values for nimbleModel()}

\item{INCLUDE_ME}{logical, default \code{TRUE}; should error in population measurements be included}

\item{seed}{}
}
\value{
list with population trajectory and observations with error if applicable; also constants and initial values used to simulate
}
\description{
Simulate population trajectories
}
