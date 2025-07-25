% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bayesian.R
\docType{class}
\name{Baysian}
\alias{Baysian}
\title{Calibrate the model to data}
\description{
A calibrator that uses a Baysian method.
}
\section{Super class}{
\code{\link[REpiSim:Calibrator]{REpiSim::Calibrator}} -> \code{Baysian}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{samples}}{the samples from the fit, a read-only field}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-Baysian-calibrate}{\code{Baysian$calibrate()}}
\item \href{#method-Baysian-clone}{\code{Baysian$clone()}}
}
}
\if{html}{\out{
<details open><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="REpiSim" data-topic="Calibrator" data-id="initialize"><a href='../../REpiSim/html/Calibrator.html#method-Calibrator-initialize'><code>REpiSim::Calibrator$initialize()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Baysian-calibrate"></a>}}
\if{latex}{\out{\hypertarget{method-Baysian-calibrate}{}}}
\subsection{Method \code{calibrate()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Baysian$calibrate(initial.values, parms, priors, guess, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{initial.values}}{the initial values for the model. The parameters 
that need to be estimate should be NA, those that do not need to be
estimated must contain a finite value.}

\item{\code{parms}}{the parameter values of the model. The parameters 
that need to be estimate should be NA, those that do not need to be
estimated must contain a finite value.}

\item{\code{priors}}{a named list of Distribution objects specifying the priors.
Each name is a parameter.}

\item{\code{guess}}{the initial guess of the parameters to be fitted}

\item{\code{...}}{extra arguments to be passed to calibrators}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Baysian-clone"></a>}}
\if{latex}{\out{\hypertarget{method-Baysian-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Baysian$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
