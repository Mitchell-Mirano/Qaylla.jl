"""
    euler(F::Function,
          a::Number,
          y_0::Number,
          b::Number,
          nh::Number)::Tuple{Vector{Float64}, Vector{Float64}}

Computes the Euler aproximation for `y(t)`, when a≤t≤b.

## Arguments
- `F(t,y)::Function`: The `y'(t)`.
- `a::Number`: The initial point.
- `y_0::Number`: The y value in a.
- `b::Number`: The final point. 
- `nh::Number`: If 1≤nh, nh is the number of subintervals beetwen a and b,
else nh is the length of the subintervals beetwen a and b.

## Return 
- `t_values`::Array{Float64}: [t₀,t₁...,tₙ], where tᵢ=a + ih
- `w_values`::Array{Float64}: [w₀,w₁...,wₙ], where wᵢ≈yᵢ=y(tᵢ)
    

## Example
```jldoctest
julia> using Qaylla

julia> F(t,y)=t-y
F (generic function with 1 method)

julia> euler(F,0,2,1,5)
([0.0, 0.2, 0.4, 0.6000000000000001, 0.8, 1.0], [2.0, 1.6, 1.32, 1.1360000000000001, 1.0288000000000002, 0.9830400000000001])

julia> euler(F,0,2,1,0.2)
([0.0, 0.2, 0.4, 0.6000000000000001, 0.8, 1.0], [2.0, 1.6, 1.32, 1.1360000000000001, 1.0288000000000002, 0.9830400000000001])
```
"""
function euler(F::Function,
               a::Number,
               y_0::Number,
               b::Number,
               nh::Number)::Tuple{Vector{Float64}, Vector{Float64}}

    if nh>=1.0
        n=nh
        h=(b-a)/n
    else
        n=(b-a)/nh
        h=nh
    end
    
    t_values::Vector{Float64}=[]
    w_values::Vector{Float64}=[]

    n::Int64=n 
    ti=a
    wi=y_0
    
    push!(t_values,ti)
    push!(w_values,wi)


    for _ in 1:n
        wi=wi+h*F(ti,wi)
        ti=ti+h
        push!(t_values,ti)
        push!(w_values,wi)
    end

    return t_values,w_values
    
end