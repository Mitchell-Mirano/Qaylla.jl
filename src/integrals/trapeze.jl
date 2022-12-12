"""
    trapeze(f::Function,a::Number,b::Number,n::Int64=1)::Float64

Computes the integral of f applying the composite trapeze method.

## Example
```jldoctest
julia> using Qaylla

julia> f(x)=sin(x)
f (generic function with 1 method)

julia> a=0
0

julia> b=pi/2
1.5707963267948966

julia> n=4
4

julia> trapeze(f,a,b,n)
0.9871158009727754
```
"""
function trapeze(f::Function,a::Number,b::Number,n::Int64=1)::Float64
    
    h=(b-a)/n
    ∑=0
    for i in 1:n-1 
        ∑+=f(a+h*i)
    end
    return (h/2)*(f(a)+2*∑+f(b))
    
end