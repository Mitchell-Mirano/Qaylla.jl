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
julia> using NumericalMethods

julia> x=[0,1,2,3,4]
5-element Vector{Int64}:
 0
 1
 2
 3
 4

julia> y= x .|> x->x^2
5-element Vector{Int64}:
 0
 1
 4
 9
 16

julia> p=lagrange(x,y)
julia> p(2.5)
6.25
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
