% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RGillespie.R
\docType{class}
\name{RGillespie}
\alias{RGillespie}
\title{R6 class implementing the Gillespie method in R
 
This is a subclass of Simulator, using an R implementation of the Gillespie method
to simulate a compartmental model.}
\description{
R6 class implementing the Gillespie method in R
 
This is a subclass of Simulator, using an R implementation of the Gillespie method
to simulate a compartmental model.

R6 class implementing the Gillespie method in R
 
This is a subclass of Simulator, using an R implementation of the Gillespie method
to simulate a compartmental model.
}
\examples{
# an SIR model
SIR = Compartmental$new(S, I, R, title="SIR")
SIR$transition(S->I ~ beta*S*I/N, N=S+I+R, name="infection")
SIR$transition(I->R ~ gamma*I, name="recovery")
g = RGillespie$new(SIR)
g$simulate(0:100, y0=c(S=1000, I=10, R=0), parms=c(beta=0.4,gamma=0.2))
}
\section{Super class}{
\code{\link[REpiSim:Simulator]{REpiSim::Simulator}} -> \code{RGillespie}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-RGillespie-clone}{\code{RGillespie$clone()}}
}
}
\if{html}{\out{
<details open><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="REpiSim" data-topic="Simulator" data-id="initialize"><a href='../../REpiSim/html/Simulator.html#method-Simulator-initialize'><code>REpiSim::Simulator$initialize()</code></a></span></li>
<li><span class="pkg-link" data-pkg="REpiSim" data-topic="Simulator" data-id="simulate"><a href='../../REpiSim/html/Simulator.html#method-Simulator-simulate'><code>REpiSim::Simulator$simulate()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-RGillespie-clone"></a>}}
\if{latex}{\out{\hypertarget{method-RGillespie-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{RGillespie$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
