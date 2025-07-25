% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Likelihood.R
\docType{class}
\name{Poisson}
\alias{Poisson}
\title{The Poisson distribution}
\description{
The Poisson distribution

The Poisson distribution
}
\details{
This distribution has no extra parameter besides the mean
}
\section{Super class}{
\code{\link[REpiSim:Likelihood]{REpiSim::Likelihood}} -> \code{Poisson}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-Poisson-logL}{\code{Poisson$logL()}}
\item \href{#method-Poisson-clone}{\code{Poisson$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Poisson-logL"></a>}}
\if{latex}{\out{\hypertarget{method-Poisson-logL}{}}}
\subsection{Method \code{logL()}}{
The log likelihood function for the Poisson distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Poisson$logL(x, mean, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{the value to calculate the likelihood}

\item{\code{mean}}{the mean of the Poisson distribution}

\item{\code{...}}{other parameters (not used)}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Poisson-clone"></a>}}
\if{latex}{\out{\hypertarget{method-Poisson-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Poisson$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
