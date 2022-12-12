"""
    trapeze2d(F::Function,
              a::Number,
              b::Number,
              c::Number,
              d::Number,
              n::Int,
              m::Int)::Float64

Computes the integral of F  in a rectangular region the plain `XY`, 
applying the composite trapeze method.

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

julia> trapeze2d(F,0,1,0,1,2,2)
0.8250805155040757

julia> trapeze2d(F,0,1,0,1,4,2)
0.7938305155040756

julia> trapeze2d(F,0,1,0,1,2,4)
0.8323009375715021

julia> trapeze2d(F,0,1,0,1,4,4)
0.8010509375715024
```
"""
function trapeze2d(F::Function,
                   a::Number,
                   b::Number,
                   c::Number,
                   d::Number,
                   n::Int,
                   m::Int)::Float64

    h=(b-a)/n
    k=(d-c)/m 

    ∫=F(a,c)+F(a,d)+F(b,c)+F(b,d)

    xᵢ::Float64=0

    for i in 1:n-1
        xᵢ=a+i*h
        ∫+=2*(F(xᵢ,c)+F(xᵢ,d))
    end

    for j in 1:m-1
        ∫+=2*(F(a,c+j*k)+F(b,c+j*k))
    end

    for j in 1:m-1
        for i in 1:n-1
            xᵢ=a+i*h
            ∫+=4*F(xᵢ,c+j*k)
        end
    end

    return (h*k/4)*∫
end

"""
    trapeze2d(F::Function,
              a::Number,
              b::Number,
              c::Function,
              d::Function,
              n::Int,
              m::Int)::Float64

Computes the integral of F  in a not rectangular region the plain `XY`, 
applying the composite trapeze method.

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

julia> trapeze2d(F,0,1,c,d,2,2)
0.7143707767687785

julia> trapeze2d(F,0,1,c,d,4,2)
0.6526600971423198

julia> trapeze2d(F,0,1,c,d,2,4)
0.7187860972222649

julia> trapeze2d(F,0,1,c,d,4,4)
0.6563466180940812
```
"""
function trapeze2d(F::Function,
                   a::Number,
                   b::Number,
                   c::Function,
                   d::Function,
                   n::Int,
                   m::Int)::Float64

    h=(b-a)/n
    k(x)=(d(x)-c(x))/m 

    ∫::Float64=k(a)*(F(a,c(a))+F(a,d(a))) + k(b)*(F(b,c(b))+F(b,d(b)))

    xᵢ::Float64=0

    for i in 1:n-1
        xᵢ=a+i*h
        ∫+=2*k(xᵢ)*(F(xᵢ,c(xᵢ))+F(xᵢ,d(xᵢ)))
    end

    for j in 1:m-1
        ∫+=2*(k(a)*F(a,c(a)+j*k(a))+k(b)*F(b,c(b)+j*k(b)))
    end

    for j in 1:m-1
        for i in 1:n-1
            xᵢ=a+i*h
            ∫+=4*k(xᵢ)*F(xᵢ,c(xᵢ)+j*k(xᵢ))
        end
    end

    return (h/4)*∫
end