% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Model.R
\docType{class}
\name{Model}
\alias{Model}
\title{R6 class representing a mathematical model.}
\description{
R6 class representing a mathematical model.

R6 class representing a mathematical model.
}
\details{
A mathematical model is described by a system of ODEs. The state variables
are called compartments, and the equations give the rate of change (i.e.,
the time derivatives of the states). An equation is specified by an
R formula, the parameters are automatically extracted. A model can then
be used to construct numerical or stochastic simulations, or be typeset
as latex equations.

The compartment name can coincide with a parameter name, 
in which case the parameter is converted into a compartment. But the
name cannot conflict with another compartment or a substitution.


If the expression uses a compartment, then the parameter is considered
a compartment when simulating the model, i.e., it is returned as a column
in the data frame returned by the simulation.


The deleted compartment or substitution is converted to a 
parameter. A parameter cannot be deleted. If the equation is changed
so that a parameter is not used by the model any more, then the paramter
is automatically removed.


The new name cannot conflict with another compartment,
substitution or paramter
}
\examples{

## ------------------------------------------------
## Method `Model$new`
## ------------------------------------------------

# An SIR model 
SIR = Model$new(
  S ~ -beta*S*I/N, # the dS/dt equation
  I ~ beta*S*I/N - gamma*I, # the dI/dt equation
  R ~ gamma*I, # the dR/dt equation
  N = S + I + R # the total population N
)
print(SIR)

## ------------------------------------------------
## Method `Model$compartment`
## ------------------------------------------------

# an SIR model
SIR = Model$new()
SIR$compartment(S ~ -beta*S*I/N)$
  compartment(I ~ beta*S*I/N - gamma*I)$
  compartment(R ~ gamma*I)$
  where(N = S + I + R)
print(SIR)

## ------------------------------------------------
## Method `Model$where`
## ------------------------------------------------

# an SIR model
SIR = Model$new()
SIR$compartment(S ~ -beta*S*I)$
  compartment(I ~ beta*S*I - gamma*I)$
  compartment(R ~ gamma*I)$
  where(pairs=list(beta=quote(b/N)))$
  where(N = S + I + R)
print(SIR)

## ------------------------------------------------
## Method `Model$delete`
## ------------------------------------------------

# an SIR model
SIR = Model$new()
SIR$compartment(S ~ -beta*S*I/N)$
  compartment(I ~ beta*S*I/N - gamma*I)$
  compartment(R ~ gamma*I)$
  where(N = S + I + R)
SIR$delete("N")
print(SIR)

## ------------------------------------------------
## Method `Model$rename`
## ------------------------------------------------

# an SIR model
SIR = Model$new()
SIR$compartment(S ~ -beta*S*I/N)$
  compartment(I ~ beta*S*I/N - gamma*I)$
  compartment(R ~ gamma*I)$
  where(N = S + I + R)
SIR$rename(S->U)$rename(beta->b)
print(SIR)
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{restricted}}{a boolean variable indicating whether to restrict functions 
allowed to be used in the ODE system. Default to FALSE}

\item{\code{attached.functions}}{the list of provided R functions to be used in
the model}
}
\if{html}{\out{</div>}}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{compartments}}{A read-only field that returns a character vector of compartment names}

\item{\code{equations}}{A read-only field that returns a named list model equations}

\item{\code{parameters}}{A read-only field that returns a character vector of parameter names}

\item{\code{substitutions}}{A read-only field that returns a named list of expressions}

\item{\code{representation}}{a read-only active field that returns the representation of the model
it returns a list that contains the equations substitutions. The compartmental model can 
then be reconstructed from the representation.}

\item{\code{order}}{a read-only active field that returns the order of the equations
and aliases that must appear to satisfy dependencies in calculation}

\item{\code{missing}}{a read-only field that returns the names of functions defined 
neither in the global environment, nor in the attached.functions list.}

\item{\code{t}}{the independent variable anme}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-Model-new}{\code{Model$new()}}
\item \href{#method-Model-compartment}{\code{Model$compartment()}}
\item \href{#method-Model-where}{\code{Model$where()}}
\item \href{#method-Model-delete}{\code{Model$delete()}}
\item \href{#method-Model-rename}{\code{Model$rename()}}
\item \href{#method-Model-format}{\code{Model$format()}}
\item \href{#method-Model-clone}{\code{Model$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Model-new"></a>}}
\if{latex}{\out{\hypertarget{method-Model-new}{}}}
\subsection{Method \code{new()}}{
constructor

It constructs a Model object with compartments and substitutions.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Model$new(..., t = "t", functions = NULL, file = NULL, .restricted = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Each extra parameter is either passed to the `compartment` method
if it is a formula with the form `name ~ value`, or passed to the
`where` methods to define a substitution if it is a named expression.}

\item{\code{t}}{the name of the independent variable, either a name or a string}

\item{\code{functions}}{a list of functions to be used in this model}

\item{\code{file}}{if not NULL, a path or connection to a model file 
to read the model from}

\item{\code{.restricted}}{a boolean variable indicating whether the ODE system
only has access to a limited set of functions. Default to TRUE}
}
\if{html}{\out{</div>}}
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{# An SIR model 
SIR = Model$new(
  S ~ -beta*S*I/N, # the dS/dt equation
  I ~ beta*S*I/N - gamma*I, # the dI/dt equation
  R ~ gamma*I, # the dR/dt equation
  N = S + I + R # the total population N
)
print(SIR)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Model-compartment"></a>}}
\if{latex}{\out{\hypertarget{method-Model-compartment}{}}}
\subsection{Method \code{compartment()}}{
define a compartment using a formula
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Model$compartment(eq)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{eq}}{a formula  in the form of `name ~ rate`. This defined a 
compartment with a given name and rate of change.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The invisible object `self` for chained operations
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{# an SIR model
SIR = Model$new()
SIR$compartment(S ~ -beta*S*I/N)$
  compartment(I ~ beta*S*I/N - gamma*I)$
  compartment(R ~ gamma*I)$
  where(N = S + I + R)
print(SIR)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Model-where"></a>}}
\if{latex}{\out{\hypertarget{method-Model-where}{}}}
\subsection{Method \code{where()}}{
define parameters as substitutions
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Model$where(..., pairs = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Each extra parameter must be named, which name is a parameter
and the value is the substitution.}

\item{\code{pairs}}{a named list, which is an alternative form to provide substitutions}
}
\if{html}{\out{</div>}}
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{# an SIR model
SIR = Model$new()
SIR$compartment(S ~ -beta*S*I)$
  compartment(I ~ beta*S*I - gamma*I)$
  compartment(R ~ gamma*I)$
  where(pairs=list(beta=quote(b/N)))$
  where(N = S + I + R)
print(SIR)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Model-delete"></a>}}
\if{latex}{\out{\hypertarget{method-Model-delete}{}}}
\subsection{Method \code{delete()}}{
delete a compartment or a substitution
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Model$delete(name)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{name}}{a character specifying the  name of the compartment 
or the substitution to be deleted.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
an invisible self to chain methods.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{# an SIR model
SIR = Model$new()
SIR$compartment(S ~ -beta*S*I/N)$
  compartment(I ~ beta*S*I/N - gamma*I)$
  compartment(R ~ gamma*I)$
  where(N = S + I + R)
SIR$delete("N")
print(SIR)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Model-rename"></a>}}
\if{latex}{\out{\hypertarget{method-Model-rename}{}}}
\subsection{Method \code{rename()}}{
rename a compartment, a substitution or a parameter
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Model$rename(formula = NULL, from = NULL, to = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{formula}}{the @param from or @param to can be specified by a formula
in the form of `from -> to` or `to <- from`.}

\item{\code{from}}{is the name of the compartment, substitution or paramter
to be changed.}

\item{\code{to}}{is the new name}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
an invisible self to chain methods
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{# an SIR model
SIR = Model$new()
SIR$compartment(S ~ -beta*S*I/N)$
  compartment(I ~ beta*S*I/N - gamma*I)$
  compartment(R ~ gamma*I)$
  where(N = S + I + R)
SIR$rename(S->U)$rename(beta->b)
print(SIR)
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Model-format"></a>}}
\if{latex}{\out{\hypertarget{method-Model-format}{}}}
\subsection{Method \code{format()}}{
format the class for printing
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Model$format()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Model-clone"></a>}}
\if{latex}{\out{\hypertarget{method-Model-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Model$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
