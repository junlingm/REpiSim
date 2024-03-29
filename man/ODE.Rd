% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ODE.R
\docType{class}
\name{ODE}
\alias{ODE}
\title{A R6 class for numerically solving the ODE models}
\description{
A R6 class for numerically solving the ODE models

A R6 class for numerically solving the ODE models
}
\details{
ODE is a subclass of `Simulator`. This class is only available if the `deSolve` 
package is installed.
}
\section{Super class}{
\code{\link[REpiSim:Simulator]{REpiSim::Simulator}} -> \code{ODE}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{model}}{a read-only field giving the R function for the ODE equations
to be passed to deSolve::ode as the `func` argument.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-ODE-new}{\code{ODE$new()}}
\item \href{#method-ODE-clone}{\code{ODE$clone()}}
}
}
\if{html}{\out{
<details open><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="REpiSim" data-topic="Simulator" data-id="simulate"><a href='../../REpiSim/html/Simulator.html#method-Simulator-simulate'><code>REpiSim::Simulator$simulate()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ODE-new"></a>}}
\if{latex}{\out{\hypertarget{method-ODE-new}{}}}
\subsection{Method \code{new()}}{
constructor
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ODE$new(model)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{model}}{the model to simulate. It must be an object of a subclass
of `Model`.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ODE-clone"></a>}}
\if{latex}{\out{\hypertarget{method-ODE-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ODE$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
