% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Likelihood.R
\name{NBinom}
\alias{NBinom}
\title{The negative binomial distribution.}
\description{
The negative binomial distribution.

The negative binomial distribution.
}
\details{
This distribution has one parameter named "size" in addition to the
mean, which is the size parameter in stats::dnbinom
}
\section{Super class}{
\code{REpiSim::Likelihood} -> \code{NBinom}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-NBinom-logL}{\code{NBinom$logL()}}
\item \href{#method-NBinom-clone}{\code{NBinom$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-NBinom-logL"></a>}}
\if{latex}{\out{\hypertarget{method-NBinom-logL}{}}}
\subsection{Method \code{logL()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NBinom$logL(x, mean, size)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-NBinom-clone"></a>}}
\if{latex}{\out{\hypertarget{method-NBinom-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NBinom$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
