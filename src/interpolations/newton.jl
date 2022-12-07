"""
    newton(x,y)::Function

Computes the Newton interpolation polynomial

## Arguments
- `x::Vector{Number}`: The xᵢ values for i=1,2...n 
- `y::Vector{Number}`: The f(xᵢ) values for i=1,2...n 

## Return 
- `p::Function`: The Newton interpolation polynomial

## Example
```jldoctest
julia> using NumericalMethods:newton

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

julia> p=newton(x,y)
julia> p(2.5)
6.25
```
"""
function newton(x,y)::Function

    n=length(x)
    N=zeros(n,n)
    
    for i in 1:n
        N[i,1]=y[i]
    end

    for i in 2:n
        for j in 2:n
            if j<=i 
                N[i,j]=(N[i,j-1]-N[i-1,j-1])/(x[i]-x[i-j+1])
            end
        end
    end
    
    function p(value::Number,progresive::Bool=true)

        p_value=0

        if progresive==true
            p_value=N[1,1]
            for i in 2:n
                pp=1
                for j in 1:i-1
                    pp*=(value-x[j])
                end
                p_value+=N[i,i]*pp
            end
        end

        if progresive==false
            p_value=N[n,1]
            for i in 2:n
                pp=1
                for j in 0:i-2
                    pp*=(value-x[n-j])
                end
                p_value+=N[n,i]*pp
            end
        end

        return p_value
    end

    return p
end