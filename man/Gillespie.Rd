% Please edit documentation in R/RGillespie.R
\name{Gillespie}
\alias{Gillespie}
\title{R6 class implementing the Gillespie method}
\description{
This is a subclass of \link{Simulator}, implementing the Gillespie stochastic
simulation algorithm for a \link{Compartmental} model. The simulator uses the
R backend by default and can use a cpp11 backend with \code{compile = TRUE}.
}
\section{Construction}{
\if{html}{\out{<div class="r">}}\preformatted{Gillespie$new(model, compile = FALSE)}\if{html}{\out{</div>}}

\describe{
\item{\code{model}}{A \link{Compartmental} model.}
\item{\code{compile}}{Logical. If \code{TRUE}, compile and use a cpp11 backend.}
}
}
\examples{
SIR <- Compartmental$new(S, I, R)
SIR$transition(S->I ~ beta*S*I/N, N=S+I+R, name="infection")
SIR$transition(I->R ~ gamma*I, name="recovery")
g <- Gillespie$new(SIR)
g$simulate(0:10, y0=c(S=1000, I=10, R=0), parms=c(beta=0.4, gamma=0.2))
}
