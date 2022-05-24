# REpiSim -- An R package for mathematical models for infectious diseases

This R package provides an interface to define compartmental infectious disease models, typeset the model in LaTeX, simulate the model numerically by solving the ODE system, ot simulate the model stochastically using the Gillespie method. If you are simply interested in simulating an ODE model, this package also provides an interface to define an ODE model with ease.

## Model Class
The R6 class R6 class to define a generic ODE model is the ```Model```. The equations may be directly defined in the ```new``` method, or individually specified using the ```compartment``` method. See more details below. The parameters of the model are automatically extracted from the equations. These parameters may also be defined as substitutions, i.e., as an expression.

The model may then be passed to a TexFormatter object to typeset as LaTeX equations, or be passed to an ODE object to numerically solve the ODE.

### Create a model

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
  parms=c(beta=0.4, gamma=0.2), 
  vars=c("I", "N")
)
```

## Compartmental class
The R6 class for defining a compartmental model is Compartmental, which is a subclass of ```Model```. The main advantage of using the Compartmental class over the Model class is that we can specify state transitions, i.e., we describe the flowchart instead of directly specifying the equations. We indeed provide a facility to generate the flowchart as the LaTeX code for a tikzpicture. In addition, the state transitions correspond to random events. We can thus run stochastic simulations using the Gillespie method.

### Create a compartmental model

The constructor of the Compartmental class takes the names of compartments. Alternatively the class names may also be specified using the ```compartment``` method, which takes one compartment name. To create the SIR model, we may start with the following statement

```
m = Compartmental$new(S, I, R)
```

The transitions are specified by the ```transition``` method using a formula
that has the form ```from -> to ~ rate``` where from and to are the compartment names, and rate is an expression for the rate of the transition. Either from or to may be NULL, specifying that new individuals is from outside the system (e.g., births) or the individuals leave the system (e.g., deaths). A transition has a name, which can be used to further modify it. It can be specified by the ```name``` argument. Alternatively, if the name is not specified, then a name is automatically generated, and is returned by the ```transition``` method. Substitutions may be passed to the transition method as named arguments. For example, the transitions for an SIR model with births and deaths may be specified as

```
# infection, where the population size N is defined as a substitution
m$transition(S->I ~ beta*S*I/N, N=S+I+R)
# recovery
m$transition(I->R ~ gamma*I)
# births
m$transition(NULL->S ~ lambda)
# deaths
death.S = m$transition(S->NULL ~ mu*S)
m$transition(I->NULL ~ mu*I)
m$transition(R->NULL ~ mu*R)
print(m)
```


The rate may be specified as a per capita rate, either by specifying ```percapita=TRUE``` in the argument, or specify the rate as ```percapita(rate)```. The following statements are equivalent.

```
m$transition(S->NULL ~ mu*S, name=death.S)
m$transition(S->NULL ~ percapita(mu), name=death.S)
m$transition(S->NULL ~ mu, percapita=TRUE, name=death.S)
```

### Stochastic simulations

Because the transitions of a compartmental model may correspond to random events that change the system state, we may use the Gillespie method to simulate the model. To do so, we use an object of the ```RGillespie``` R6 class, which is an R implementation. Like the ODE simulate, the constructor takes a compartmental model (note that it does not work with the Model class). To simulate, we use the ```simulate``` method, like the ODE simulation.

```
sim = RGillespie$new(m)
x = sim$simulate(
  t=0:100, 
  y0=c(S=1000, I=10, R=0),
  parms=c(beta=0.4, gamma=0.2, lambda=100, mu=0.1))
print(head(x))
```
