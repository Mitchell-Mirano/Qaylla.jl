"""
    simpson(f::Function,a::Number,b::Number,n::Int64=2)::Float64

Computes the integral of f applying the composite Simpson method.

## Example
```jldoctest
julia> using NumericalMethods

julia> f(x)=sin(x)
f (generic function with 1 method)

julia> a=0
0

julia> b=pi/2
1.5707963267948966

julia> n=4
4

julia> simpson(f,a,b,n)
1.0001345849741938
```
"""
function simpson(f::Function,a::Number,b::Number,n::Int64=2)::Float64
    
    if n%2!=0
        n=n+1
    end
    h=(b-a)/n
    m=Integer(n/2)

    ∑1=0
    for i in 1:m
        ∑1+=f(a+(2*i-1)*h)
    end

    ∑2=0
    for i in 1:m-1
        ∑2+=f(a+(2*i)*h)
    end

    return (h/3)*(f(a)+f(b)+4*∑1+2*∑2)
    
end