"""
    lagrange(x,y)::Function

Computes the Lagrange interpolation polynomial

## Arguments
- `x::Vector{Number}`: The xᵢ values for i=1,2...n 
- `y::Vector{Number}`: The f(xᵢ) values for i=1,2...n 

## Return 
- `p::Function`: The Lagrange interpolation polynomial

## Example
```jldoctest
using Yaqa

x=0:0.5:3
y=exp.(x)
p=newton(x,y)
p(2.5)

# output

12.182493960703471
```
"""
function lagrange(x,y)::Function
    
    function p(value::Number)
        p_value=0
        n=length(x)

        for i in 1:n
            li=1
            for j in 1:n
                if i!=j
                    li*=(value-x[j])/(x[i]-x[j])
                end
            end
            p_value+=li*y[i]
        end
        
        return p_value
    end

    return p
end 
