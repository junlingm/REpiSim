# REpiSim — An R package for mathematical models for infectious diseases

This R package provides an interface to define compartmental infectious disease models, typeset the model in LaTeX, simulate the model numerically by solving systems of ODEs, or simulate the model stochastically using the Gillespie method. If you are simply interested in simulating an ODE model, this package also provides a convenient interface to define an ODE model with ease.

REpiSim supports two types of models:

- Generic **ODE models**, defined directly by differential equations.
- **Compartmental models**, defined by transitions between states.

---

## Model class

The R6 class used to define a generic ODE model is `Model`.  
Equations may be defined either:

- directly in the `new()` method, or
- added one at a time using the `compartment()` method.

Parameters appearing in the equations are automatically detected. Algebraic substitutions can be defined using the `where()` method.

### Creating a model

An equation is specified using R formula syntax. For example,

```r
S ~ -beta*S*I
```

defines a compartment `S` with derivative

\[
\frac{dS}{dt} = -\beta S I.
\]

A basic SIR model can be defined as follows:

```r
m <- Model$new(
  S ~ -beta*S*I/N,
  I ~  beta*S*I/N - gamma*I,
  R ~  gamma*I
)
print(m)
```

Here `N` appears as a parameter. It can instead be defined as a substitution:

```r
m$where(N = S + I + R)
print(m)
```

---

## Typesetting equations in LaTeX

Differential equations and substitutions can be typeset using a `TexFormatter` object:

```r
tex <- TexFormatter$new()
tex$typeset.equation(m$equations$equations)
tex$typeset.equation(m$equations$where)
```

Individual expressions can also be typeset:

```r
tex$typeset(quote(beta*S*I))
```

Custom LaTeX symbols can be defined using `set.symbols()`.

---

## Numerical simulation

Numerical simulation of an ODE model is handled by the `ODE` class, which wraps the `deSolve` package.

```r
sim <- ODE$new(m)
x <- sim$simulate(
  t     = 0:100,
  y0    = c(S = 1000, I = 10, R = 0),
  parms = c(beta = 0.4, gamma = 0.2)
)
head(x)
```

The returned object is a data frame containing the simulated trajectories. You may select a subset of variables via the `vars` argument, and substitution values can also be returned.

---

## Compartmental models

For many infectious disease models, it is more natural to describe the dynamics in terms of state transitions (for example, Susceptible → Infected → Recovered) rather than writing differential equations explicitly.

REpiSim provides the `Compartmental` class for this purpose. The `Compartmental` class is a subclass of `Model`. Internally, all transition rules are automatically translated into a system of ODEs, so everything that works for `Model` (LaTeX typesetting, numerical simulation, etc.) also works for `Compartmental`.

### Defining compartments

Compartments (state variables) are declared when creating the model:

```r
m <- Compartmental$new(S, I, R)
```

This declares `S`, `I`, and `R` as state variables, with their derivatives initially set to zero.

### Adding transitions

Transitions are added using the `transition()` method. A transition has the general form:

```r
from -> to ~ rate
```

A standard SIR model can be written as:

```r
m$transition(S -> I ~ beta*S*I/N, N = S + I + R)
m$transition(I -> R ~ gamma*I)
```

This corresponds to the ODE system

\[
\begin{aligned}
\frac{dS}{dt} &= -\beta S I / N, \\
\frac{dI}{dt} &= \beta S I / N - \gamma I, \\
\frac{dR}{dt} &= \gamma I.
\end{aligned}
\]

Model parameters such as `beta` and `gamma` are automatically detected. Substitutions (for example `N = S + I + R`) may be supplied directly in the transition call.

### Births, deaths, and external flows

Births, deaths, or other external flows can be modeled using `NULL` as the source or destination:

- Births into a compartment:
  ```r
  m$transition(NULL -> S ~ mu)
  ```
- Deaths from a compartment:
  ```r
  m$transition(I -> NULL ~ mu*I)
  ```

### Per-capita transition rates

For per-capita transition rates, you may either write them explicitly or use the `percapita` option:

```r
m$transition(I -> R ~ gamma, percapita = TRUE)
```

This is interpreted as a transition rate equal to `gamma * I`.

### Inspecting and simulating compartmental models

Because `Compartmental` is a subclass of `Model`, the generated differential equations can be inspected directly:

```r
print(m)
```

They can also be typeset using `TexFormatter`.

Once defined, compartmental models can be simulated in the same way as generic models.

Deterministic (ODE) simulation:

```r
sim <- ODE$new(m)
x <- sim$simulate(
  t     = 0:100,
  y0    = c(S = 1000, I = 10, R = 0),
  parms = c(beta = 0.4, gamma = 0.2)
)
```

Stochastic (Gillespie) simulation:

```r
sim <- RGillespie$new(m)
x <- sim$simulate(
  t     = 0:100,
  y0    = c(S = 1000, I = 10, R = 0),
  parms = c(beta = 0.4, gamma = 0.2)
)
```

### Choosing between `Model` and `Compartmental`

- Use **`Model`** when you want full control over the differential equations or are working with non-compartmental systems.
- Use **`Compartmental`** when your model is naturally described by transitions between states; this approach is often more concise and less error-prone.

---

## Flowcharts

REpiSim can generate TikZ flowcharts for compartmental models using the `flowchart()` function.

---

## Summary

REpiSim provides tools to:

- define ODE models using symbolic equations,
- define compartmental models using transition rules,
- typeset equations in LaTeX,
- simulate models deterministically and stochastically, and
- generate flowcharts for visualization.
