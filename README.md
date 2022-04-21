# REpiSim -- An R package for mathematical models for infectious diseases

This R package provides an interface to define compartmental infectious disease models, typeset the model in LaTeX, simulate the model numerically by solving the ODE system, ot simulate the model stochastically using the Gillespie method. If you are simply interested in simulating an ODE model, this package also provides an interface to define an ODE model with ease.

## Compartmental Class
The main interface for defining a compartmental model is the R6 class named Compartmental. The constructor takes compartmental names as input, and a title (the name of the model). Then you can specify the transitions (flows) between compartments using the transition method. The parameters of the model are automatically extracted from the transition rates. These parameters may also be defined as substitutions, i.e., as an expression.

The model may then be passed to a TexFormatter object to typeset as LaTeX equations, or be passed to an ODE object to numerically solve the ODE, or be passed to an RGillespie object for stochastic simulation (you may also use CGillespie to do the simulation using C++, which is about 10x faster).
