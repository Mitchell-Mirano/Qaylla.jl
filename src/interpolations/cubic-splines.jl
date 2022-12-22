"""
    cubic_splines(x_data,y_data)::Function

Computes the cubic splines  interpolation polynomial

## Arguments
- `x_data`: The xᵢ values for i=1,2...n 
- `y_data`: The f(xᵢ) values for i=1,2...n 

## Return 
- `p::Function`: The cubic splines interpolation polynomial

## Example
```jldoctest
using Qaylla

x=0:pi/4:3*pi
y=sin.(x)
p=cubic_splines(x,y)
p(2.5)

# output

0.59842733419271

```
"""
function cubic_splines(x_data,y_data)::Function
    n=length(x_data)
    hi::Array{Float64}=zeros(n-1)
    A::Matrix{Float64}=zeros(n,n)
    b::Array{Float64}=zeros(n)
    ai::Array{Float64}=y_data[1:n-1]
    bi::Array{Float64}=zeros(n)
    ci::Array{Float64}=zeros(n)
    di::Array{Float64}=zeros(n)

    splines::Array{Function}=[]
    
    for i in 1:n-1
        hi[i]=x_data[i+1]-x_data[i]
    end

    A[1,1]=1
    A[n,n]=1

    for i in 2:n-1
        A[i,i-1]=hi[i-1]
        A[i,i]=2*(hi[i-1]+hi[i])
        A[i,i+1]=hi[i]
    end

    b[1]=0
    b[n]=0
    for i in 2:n-1
        b[i]=(3/hi[i])*(y_data[i+1]-y_data[i])-(3/hi[i-1])*(y_data[i]-y_data[i-1])
    end

    ci=inv(A)*b

    for i in 1:n-1
        bi[i]=(1/hi[i])*(y_data[i+1]-y_data[i])-(hi[i]/3)*(2*ci[i]+ci[i+1])
    end 

    for i in 1:n-1
        di[i]=(ci[i+1]-ci[i])/(3*hi[i])
    end

    for i in 1:n-1
        si(x)=ai[i] + bi[i]*(x-x_data[i]) + ci[i]*(x-x_data[i])^2 + di[i]*(x-x_data[i])^3
        push!(splines,si)
    end

    function p(value::Number)::Float64
        
        value_interval::Int64=1
        
        for j in 1:n-1
            if x_data[j]<=value && value<=x_data[j+1]
                value_interval=j
                break
            end 
        end
        return splines[value_interval](value)
    end

    return p
end 

