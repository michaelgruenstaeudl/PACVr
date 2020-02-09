\name{RCircos.Env}
\alias{RCircos.Env}

\title{
Export the custom environment 'RCircos.Env'
}

\description{
This function exports the custom environment 'RCircos.Env' used by \pkg{RCircos}.
}

\details{
\pkg{PACVr} employs \pkg{RCircos} as its visualization engine. In its operation, \pkg{RCircos} defines a custom environment (called 'RCircos.Env') and reads/writes variables to this environment from various of its functions. In order to make this environment accessible to \pkg{RCircos} within \pkg{PACVr} and, simultaneously, fulfil the requirements of CRAN, this export command was created. For more information, please see the Stackoverflow post at \url{https://stackoverflow.com/questions/56875962/r-package-transferring-environment-from-imported-package}.
}
