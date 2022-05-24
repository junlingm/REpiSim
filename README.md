# REpiSim -- An R package for mathematical models for infectious diseases

This R package provides an interface to define compartmental infectious disease models, typeset the model in LaTeX, simulate the model numerically by solving the ODE system, ot simulate the model stochastically using the Gillespie method. If you are simply interested in simulating an ODE model, this package also provides an interface to define an ODE model with ease.

## Model Class
The main interface to define a generic ODE model is the ```Model``` R6 class. The equations may be directly defined in the ```new``` method, or individually specified using the ```compartment``` method. See more details below. The parameters of the model are automatically extracted from the equations. These parameters may also be defined as substitutions, i.e., as an expression.

The model may then be passed to a TexFormatter object to typeset as LaTeX equations, or be passed to an ODE object to numerically solve the ODE.

### Creating a model

To create a model, we can use the ```Model$new()``` method with either no argument (which createa an empty model) and then use the ```compartment``` method to define each class, or directly pass the equations to the ```new``` method.

An equation can be expressed as a formula. For example, `S ~ - beta*S*I` defines a compartment S with the equation
\[
\frac{dS}{dt} = \beta S I
\]
The following code defines a SIR model with standard incidence.
```
m = Model$new(
  S ~ -beta*S*I/N,
  I ~ beta*S*I/N - gamma*I,
  R ~ gamma*I
)
print(m)
```
It shows that there are three classes, S, I and R, and it uses three parameters, beta, N and gamma.

Here N is a free parameter, but it should be the total population. We can define a parameter as a substitution, using the ```where``` method, using an argument with the name N and the expression ```S+I+R``` as the value, i.e.,

```
m$where(N=S+I+R)
print(m)
```


### Typeset the equations in LaTeX
A model can be typeset as LaTeX equations using an object of the ```TexFormatter``` class. A TexFormatter obect has a typeset.equation method that which accepts a list of R expressions. For example, the ```equations```
field of a Model object contains two of such lists, the ```equations``` component contains the differential equations, and the ```where``` component contains the substitutions.

```
tex = TexFormatter$new()
tex$typeset.equation(m$equations$equations)
tex$typeset.equation(m$equations$where)
```

A TexFormatter object also has a typeset method which typesets an R expression as LaTeX.

```
tex$typeset(quote(beta*S*I))
```

A TexFormatter object automatically recognize greek letters trig functions and \(\infty\). We may define new symbols using the ```set.symbols``` method, which takes a named list with the names as the R names and the values latex commands. For example, the population size N in the SIR model above is not a compartment. To make it clear, we will use the corresponding latex mathcal symbol.

```
tex$set.symbols(list(N="\\mathcal{N}"))
tex$typeset.equation(m$equations$equations)
```

### Numerical simulations

To numerically simulate the ODE model, we need an object of class ```ODE```. Note that we need the ```deSolve``` package installed in order to use the class. The constructor of the ODE class takes the Model object to simulate. The actual simulation is done using the ```simulate``` class, which takes a numeric vector of time points to simulate, the initial conditions, and the parameter values. Additional arguments are passed unchanged to the underlying ```deSolve::ode``` function. The simulate method returns a data.frame which first column is the time, and the following columns are the state values. The following example simulates the SIR model defined above.

```
sim = ODE$new(m)
x = sim$simulate(
  t=0:100, 
  y0=c(S=1000, I=10, R=0), 
  parms=c(beta=0.4, gamma=0.2))
head(x)
```

We may also select a subset of state variables to simulate, using the ```vars``` argument of the ```simulate``` method. It takes a vector of characters giving the variable names. We may also ask it to return the value of a substitution. The following example returns both I and N of the SIR model so that we can verify that the total population N is indeed a constant.

```
sim$simulate(
  t=0:100, 
  y0=c(S=1000, I=10, R=0), 
  parms=c(beta=0.4, gamma=0.2), vars=c("I", "N")
)
```

## Compartmental class
The main interface for defining a compartmental model is the R6 class named Compartmental, which is a subclass of ```Model```. The constructor takes compartmental names as input. Then you can specify the transitions (flows) between compartments using the transition method. The parameters of the model are automatically extracted from the transition rates. These parameters may also be defined as substitutions, i.e., as an expression.

Because it is a subclass of Model, we may use TexFormatter and ODE with it. In addition, a Compartmental object may be passed to an RGillespie object for stochastic simulation. You may also use CGillespie to do the simulation using C++, which is about 10x faster but slower in compilation. The CGillespie class needs the ```cpp11``` package to be installed.

