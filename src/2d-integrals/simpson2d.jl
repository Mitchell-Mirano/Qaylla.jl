"""
    simpson2d(F::Function,
              a::Number,
              b::Number,
              c::Number,
              d::Number,
              n::Int,
              m::Int)::Float64

Computes the integral of F  in a rectangular region the plain `XY`, 
applying the composite Simpson method.

## Arguments
- `F(x,y)::Function`: F(x,y) in the function to integrate.
- `a::Number`: The lower limit of the region in x.
- `b::Number`: The upper limit of the region in x.
- `c::Number`: The lower limit of the region in y.
- `d::Number`: The upper limit of the region in y.
- `n::Number`: The number of subintervals beetwen a and b.
- `m::Number`: The number of subintervals beetwen c and d.


## Example
```jldoctest
julia> using Yaqa

julia> F(x,y)=x^2 + sin(y)
F (generic function with 1 method)

julia> simpson2d(F,0,1,0,1,2,2)
0.793195523204118

julia> simpson2d(F,0,1,0,1,2,4)
0.7930410782606442

julia> simpson2d(F,0,1,0,1,4,2)
0.793195523204118

julia> simpson2d(F,0,1,0,1,4,4)
0.7930410782606443

julia> simpson2d(F,0,1,0,1,100,100)
0.7930310274907315
```
"""
function simpson2d(F::Function,
                   a::Number,
                   b::Number,
                   c::Number,
                   d::Number,
                   n::Int,
                   m::Int)::Float64

    if n%2!=0
        n+=1
    end

    if m%2!=0
        m+=1
    end

    h=(b-a)/n
    k=(d-c)/m 
    n_mean::Int=n/2
    m_mean::Int=m/2

    ∫=F(a,c)+F(a,d)+F(b,c)+F(b,d)

    for i in 1:n_mean-1
        ∫+=2*(F(a+h*2*i,c) + F(a+h*2*i,d))
    end

    for j in 1:m_mean-1
        ∫+=2*(F(a,c+k*2*j) + F(b,c+k*2*j))
    end 

    for i in 1:n_mean
        ∫+=4*(F(a+h*(2*i-1),c) + F(a+h*(2*i-1),d))
    end 

    for j in 1:m_mean
        ∫+=4*(F(a,c+k*(2*j-1)) + F(b,c+k*(2*j-1)))
    end

    for j in 1:m_mean-1
        for i in 1:n_mean-1
            ∫+=4*F(a+h*2*i,c+k*2*j)
        end
    end

    for j in 1:m_mean-1
        for i in 1:n_mean
            ∫+=8*F(a+h*(2*i-1),c+k*2*j)
        end
    end 

    for j in 1:m_mean
        for i in 1:n_mean-1
            ∫+=8*F(a+h*2*i,c+k*(2*j-1))
        end 
    end

    for j in 1:m_mean
        for i in 1:n_mean
            ∫+=16*F(a+h*(2*i-1),c+k*(2*j-1))
        end 
    end 

    return (h*k/9)*∫
    
end

"""
    simpson2d(F::Function,
              a::Number,
              b::Number,
              c::Function,
              d::Function,
              n::Int,
              m::Int)::Float64

Computes the integral of F  in a not rectangular region the plain `XY`, 
applying the composite Simpson method.

## Arguments
- `F(x,y)::Function`: F(x,y) in the function to integrate.
- `a::Number`: The lower limit of the region in x.
- `b::Number`: The upper limit of the region in x.
- `c(x)::Function`: The lower limit of the region in y.
- `d(x)::Function`: The upper limit of the region in y.
- `n::Number`: The number of subintervals beetwen a and b.
- `m::Number`: The number of subintervals beetwen c(x) and d(x).


## Example
```jldoctest
julia> using Yaqa

julia> F(x,y)=x^2 + sin(y)
F (generic function with 1 method)

julia> c(x)=x
c (generic function with 1 method)

julia> d(x)=2*x
d (generic function with 1 method)

julia> simpson2d(F,0,1,c,d,2,2)
0.6343236523627963

julia> simpson2d(F,0,1,c,d,4,2)
0.6367323875670757

julia> simpson2d(F,0,1,c,d,2,4)
0.6342654852512393

julia> simpson2d(F,0,1,c,d,4,4)
0.6366813209795263
```
"""
function simpson2d(F::Function,
                   a::Number,
                   b::Number,
                   c::Function,
                   d::Function,
                   n::Int,
                   m::Int)::Float64

    if n%2!=0
        n+=1
    end

    if m%2!=0
        m+=1
    end

    h::Float64=(b-a)/n
    k(x)=(d(x)-c(x))/m 

    n_mean::Int=n/2
    m_mean::Int=m/2

    ∫::Float64=k(a)*(F(a,c(a))+F(a,d(a)))+ k(b)*(F(b,c(b))+F(b,d(b)))

    xᵢ::Float64=0.0

    for i in 1:n_mean-1
        xᵢ=a+h*2*i
        ∫+=2*k(xᵢ)*(F(xᵢ,c(xᵢ)) + F(xᵢ,d(xᵢ)))
    end


    for j in 1:m_mean-1
        ∫+=2*(k(a)*F(a,c(a)+k(a)*2*j) + k(b)*F(b,c(b)+k(b)*2*j))
    end 

    for i in 1:n_mean
        xᵢ=a+h*(2*i-1)
        ∫+=4*k(xᵢ)*(F(xᵢ,c(xᵢ)) + F(xᵢ,d(xᵢ)))
    end 

    for j in 1:m_mean
        ∫+=4*(k(a)*F(a,c(a)+k(a)*(2*j-1)) + k(b)*F(b,c(b)+k(b)*(2*j-1)))
    end

    for j in 1:m_mean-1
        for i in 1:n_mean-1
            xᵢ=a+h*2*i
            ∫+=4*k(xᵢ)*F(xᵢ,c(xᵢ)+k(xᵢ)*2*j)
        end
    end

    for j in 1:m_mean-1
        for i in 1:n_mean
            xᵢ=a+h*(2*i-1)
            ∫+=8*k(xᵢ)*F(xᵢ,c(xᵢ)+k(xᵢ)*2*j)
        end
    end 

    for j in 1:m_mean
        for i in 1:n_mean-1
            xᵢ=a+h*2*i
            ∫+=8*k(xᵢ)*F(xᵢ,c(xᵢ)+k(xᵢ)*(2*j-1))
        end 
    end

    for j in 1:m_mean
        for i in 1:n_mean
            xᵢ=a+h*(2*i-1)
            ∫+=16*k(xᵢ)*F(xᵢ,c(xᵢ)+k(xᵢ)*(2*j-1))
        end 
    end 

    return (h/9)*∫

end