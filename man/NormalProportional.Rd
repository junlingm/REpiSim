% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Likelihood.R
\name{NormalProportional}
\alias{NormalProportional}
\title{A normal distribution which standard deviation is proportion to mean}
\description{
A normal distribution which standard deviation is proportion to mean

A normal distribution which standard deviation is proportion to mean
}
\details{
This distribution has one parameter named "ceof" in addition to the
mean, where sd = coef*mean
}
\section{Super class}{
\code{REpiSim::Likelihood} -> \code{NormalProportional}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-NormalProportional-logL}{\code{NormalProportional$logL()}}
\item \href{#method-NormalProportional-clone}{\code{NormalProportional$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-NormalProportional-logL"></a>}}
\if{latex}{\out{\hypertarget{method-NormalProportional-logL}{}}}
\subsection{Method \code{logL()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NormalProportional$logL(x, mean, coef)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-NormalProportional-clone"></a>}}
\if{latex}{\out{\hypertarget{method-NormalProportional-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NormalProportional$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
