# Integrales

## trapeze method
```@docs
trapeze(f::Function,a::Number,b::Number,n::Int64=1)
```
## Simpson method
```@docs
simpson(f::Function,a::Number,b::Number,n::Int64=2)
```
## Romberg method
```@docs
romberg(f::Function,a::Number,b::Number,n::Int64=4)
```
## Newton-Cotes methods
```@docs
newton_cotes(f::Function,a::Number,b::Number,n::Int64,closed::Bool=true)
```
## trapeze2d method
```@docs
trapeze2d(F::Function,a::Number,b::Number,c::Number,d::Number,n::Int,m::Int)
```
## Simpson2d method
```@docs
simpson2d(F::Function,a::Number,b::Number,c::Number,d::Number,n::Int,m::Int)
```