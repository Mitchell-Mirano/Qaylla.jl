# NumericalMethods.jl

Interpolate, derive, integrate and solve differential equations, using the main methods of numerical analysis, with the Julia power.


### Author
- [Mitchell Mirano Caro](https://github.com/Mitchell-Mirano), Faculty of Mathematics, National University of San Marcos(UNMSM).

### License
The NumericalMethods.jl package is licensed under the MIT  License.

# Installation
```julia
using Pkg
Pkg.add("https://github.com/Mitchell-Mirano/NumericalMethods.jl.git")
```
# Examples

## Interpolations

```julia
x=0:0.5:3
y=x .|> x->exp(x)
p=lagrange(x,y) #computes p(x), the Lagrange polynomial interpolation.
p(2.1) # use p(x) in a point.
```
```julia
x=0:0.5:3
y=x .|> x->exp(x)
p=newton(x,y) #computes p(x), the Newton polynomial interpolation.
p(2.1) # use p(x) in a point.
```

## Integrations

Use Simpson's composite method to compute the integral of f(x), when x is between a and b, for example:

$$\int_{0}^{\pi/2}\sin(x)dx$$

```julia
f(x)=sin(x)
a=0
b=pi/2
n=4
simpson(f,a,b,n)
```
Use Simpson's composite method to compute the integral of F(x,y), in an irregular region, when x is between a and b and y is between c(x) and d(x), for example:

$$\int_{0}^{1}\int_{x}^{2x}x^{2} + \sin(y) dydx$$

```julia
F(x,y)=x^2 + sin(y)
a,b=0,1
c(x)=x
d(x)=2*x
n,m=4,4
simpson2d(F,a,b,c,d,n,m)
```
# Development Status
This package is currently under development and will soon add new features.

If you want to participate in the maintenance or active development of the package, feel free to get in touch via an issue on GitHub or by writing an email to the developer.

# Contributing
All help is always welcome. You can help by improving the documentation, reporting bugs or adding new methods.
