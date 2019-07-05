\name{RCircos.Env}
\alias{RCircos.Env}

\title{
Exporting the custom environment used by RCircos
}

\description{
This function exports the custom environment used by RCircos.
}

\details{
One of the peculiarities of RCircos is that it defines a custom environment (called RCircos.Env) and reads/writes variables to this environment from various of its functions. For example, function RCircos.Initialize.Plot.Parameters reads and writes to this environment. In order to make this enviornment accessible to RCircos within PACVr and simultaneously fulfil the requirements of CRAN, this export command was created. Please see the following Stackoverflow post for details: \url{https://stackoverflow.com/questions/56875962/r-package-transferring-environment-from-imported-package}
}
