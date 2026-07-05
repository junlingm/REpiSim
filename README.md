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

Compartments (state variables) may be declared when creating the model:

```r
m <- Compartmental$new(S, I, R)
```

This declares `S`, `I`, and `R` as state variables, with their derivatives initially set to zero. Explicit declarations are useful when you want to control output order or include a compartment before adding transitions.

Compartment declarations are optional. If a transition mentions a compartment that has not yet been declared, it is added automatically:

```r
m <- Compartmental$new()
m$transition(I <- S ~ beta*S*I/N, N = S + I + R)
m$transition(R <- I ~ gamma, percapita = TRUE)
```

### Adding transitions

Transitions are added using the `transition()` method. A transition has the general form:

```r
to <- from ~ rate
```

A standard SIR model can be written as:

```r
m$transition(I <- S ~ beta*S*I/N, N = S + I + R)
m$transition(R <- I ~ gamma*I)
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
  m$transition(S <- NULL ~ mu)
  ```
- Deaths from a compartment:
  ```r
  m$transition(NULL <- I ~ mu*I)
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
sim <- Gillespie$new(m)
x <- sim$simulate(
  t     = 0:100,
  y0    = c(S = 1000, I = 10, R = 0),
  parms = c(beta = 0.4, gamma = 0.2)
)
```

Use `Gillespie$new(m, compile = TRUE)` to run the same simulator through the
optional cpp11 backend.

### Choosing between `Model` and `Compartmental`

- Use **`Model`** when you want full control over the differential equations or are working with non-compartmental systems.
- Use **`Compartmental`** when your model is naturally described by transitions between states; this approach is often more concise and less error-prone.

---

## Stratified models

REpiSim supports stratified state variables in both `Model` and `Compartmental`. A stratified compartment is written using indexed notation such as `S[i]`, where the index values are supplied by an index set.

For example, a two-group SIR model can be written directly as ODE equations:

```r
groups <- c("A", "K")

m <- Model$new(
  S[i] ~ -Sum(beta[i, j] * I[j], j = groups) * S[i],
  I[i] ~  Sum(beta[i, j] * I[j], j = groups) * S[i] - gamma[i] * I[i],
  R[i] ~  gamma[i] * I[i],
  .index = list(i = groups, j = groups)
)
```

The model-level equations remain readable in indexed form, while simulators receive a flat internal representation with variables such as `S_A`, `S_K`, `I_A`, and `I_K`.

The same model can be written as stratified transitions:

```r
m <- Compartmental$new(.index = list(i = groups, j = groups))

m$transition(
  I[i] <- S[i] ~ Sum(beta[i, j] * I[j], j = groups),
  percapita = TRUE,
  name = "infection"
)

m$transition(
  R[i] <- I[i] ~ gamma[i],
  percapita = TRUE,
  name = "recovery"
)
```

This creates separate flat event transitions, such as `infection_A`, `infection_K`, `recovery_A`, and `recovery_K`. That preserves the event structure needed by `Gillespie`.

### Initial values and parameters

Initial values may be supplied either in flat form:

```r
y0 <- c(
  S_A = 0.5, S_K = 0.5,
  I_A = 1e-4, I_K = 1e-4,
  R_A = 0, R_K = 0
)
```

or as a structured list:

```r
y0 <- list(
  S = c(A = 0.5, K = 0.5),
  I = c(A = 1e-4, K = 1e-4),
  R = c(A = 0, K = 0)
)
```

Parameters may also be structured:

```r
beta <- matrix(
  c(2.0, 1.5,
    1.5, 2.0),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(groups, groups)
)

parms <- list(
  beta = beta,
  gamma = c(A = 0.2, K = 0.2)
)
```

Then simulate as usual:

```r
sim <- ODE$new(m)
x <- sim$simulate(
  t = seq(0, 50, by = 0.1),
  y0 = y0,
  parms = parms
)
```

The simulation output is always a data frame with flat variable columns.

For performance, simulator functions use integer indexing into structured parameters. When names or dimnames are present, REpiSim checks that they are in the model's stratum order before simulation.

### Calibration with stratified models

Calibration accepts the same flat and structured forms for `initial.values`, `parms`, and `guess`.

For example, fixed and fitted values can be mixed using `NA` placeholders:

```r
fit$calibrate(
  initial.values = list(
    S = c(A = 0.5, K = 0.5),
    I = c(A = NA, K = 1e-4)
  ),
  parms = list(
    beta = matrix(c(1.5, NA, NA, 2.0), nrow = 2),
    gamma = c(A = 0.2, K = 0.2)
  ),
  guess = list(
    I = c(A = 1e-4),
    beta = matrix(c(NA, 1.6, 1.7, NA), nrow = 2)
  )
)
```

Equivalently, flat component names can be used:

```r
fit$calibrate(
  initial.values = c(S_A = 0.5, S_K = 0.5, I_A = 1e-4, I_K = 1e-4),
  parms = list(beta_A_A = 1.5, beta_K_K = 2.0, gamma_A = 0.2, gamma_K = 0.2),
  guess = c(beta_K_A = 1.6, beta_A_K = 1.7)
)
```

---

## Calibration

REpiSim provides a common calibration interface through the `Calibrator` class and its subclasses. A calibrator combines:

- a model,
- observed data,
- a simulator,
- a mapping from data columns to model variables, and
- a rule for estimating unknown initial values and parameters.

The main user-facing calibrators are:

- `MLE`, for maximum likelihood estimation using `bbmle`,
- `Metrop`, for Bayesian calibration using `mcmc::metrop`, and
- `LeastSquare`, for least-squares fitting.

### Likelihood-based calibration

As a small example, suppose an observed quantity decays exponentially:

```r
m <- Model$new(S ~ -a*S)

dat <- data.frame(
  time = 0:5,
  obs  = c(10.0, 8.7, 7.6, 6.5, 5.8, 5.0)
)
```

The data column `obs` is mapped to the model variable `S` using a formula:

```r
fit <- MLE$new(
  model      = m,
  time       = "time",
  data       = dat,
  likelihood = Normal(sd = 1),
  obs ~ S,
  CI = FALSE
)

fit$calibrate(
  initial.values = c(S = 10),
  parms          = numeric(),
  guess          = c(a = 0.1)
)
```

In `initial.values` and `parms`, numeric values are fixed. Omitted values are fitted. Thus, in the example above, `S(0)` is fixed at 10 and `a` is estimated.

Set `compile = TRUE` in `MLE`, `Metrop`, or `LeastSquare` to use the compiled
ODE backend during calibration.

### Fitting transformed parameters

Parameters may also be supplied as quoted expressions. This is useful for enforcing constraints. For example, to estimate an unconstrained parameter `b` while forcing `a > 0`, write:

```r
fit$calibrate(
  initial.values = c(S = 10),
  parms          = list(a = quote(exp(b))),
  guess          = c(b = log(0.1))
)
```

The calibrator fits `b`, then computes `a = exp(b)` before simulating the model.

The same pattern works for likelihood parameters. For example, a normal likelihood can estimate `log_sd` while using `sd = exp(log_sd)`:

```r
fit <- MLE$new(
  model      = m,
  time       = "time",
  data       = dat,
  likelihood = Normal(sd = quote(exp(log_sd))),
  obs ~ S,
  CI = FALSE
)

fit$calibrate(
  initial.values = c(S = 10),
  parms          = numeric(),
  guess          = c(a = 0.1, log_sd = 0)
)
```

### Distributions

Distributions are declared by their canonical parameters and density functions. The built-in distributions can be used either as fully specified prior distributions or as mean-parameterized likelihoods.

For example:

```r
Normal(mean = 0, sd = 1)  # fully specified density
Normal(sd = 1)            # likelihood: mean comes from the model
Poisson()                 # likelihood: lambda is the model mean
NBinom(size = 10)         # likelihood: prob is computed from mean and size
```

### Bayesian calibration

`Metrop` uses the same data/model/mapping interface as `MLE`, but calibration additionally requires priors and an initial guess.

```r
bfit <- Metrop$new(
  model      = m,
  time       = "time",
  data       = dat,
  likelihood = Normal(sd = quote(exp(log_sd))),
  obs ~ S
)

bfit$calibrate(
  initial.values = c(S = 10),
  parms          = list(a = quote(exp(b))),
  priors         = list(
    b      = Normal(mean = log(0.1), sd = 1),
    log_sd = Normal(mean = 0, sd = 1)
  ),
  guess = c(b = log(0.1), log_sd = 0),
  nbatch = 1000
)
```

Prior distributions can depend on previously defined parameters by using quoted expressions. For example, the following prior specification enforces `a < b`:

```r
priors <- list(
  a = Uniform(min = 0, max = 1),
  b = Uniform(min = quote(a), max = 1)
)
```

When evaluating the prior for `b`, the calibrator first uses the current value of `a` to compute the lower bound.

---

## Flowcharts

REpiSim can generate TikZ flowcharts for compartmental models using the `flowchart()` function.

---

## Summary

REpiSim provides tools to:

- define ODE models using symbolic equations,
- define compartmental models using transition rules,
- typeset equations in LaTeX,
- simulate models deterministically and stochastically,
- calibrate models using likelihood-based or Bayesian methods, and
- generate flowcharts for visualization.
