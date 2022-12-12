"""
    euler(F::Function,
          t₀::Number,
          y₀::Number,
          tₙ::Number,
          nh::Number,
          return_values::Bool=false)

Computes the Euler aproximation for `y(t)` in `tₙ`.

## Arguments
- `F(t,y)::Function`: The `y'(t)`.
- `t₀::Number`: The initial point. 
- `tₙ::Number`: The final point. 
- `nh::Number`: If 1≤nh, nh is the number of subintervals beetwen t₀ and tₙ,
else nh is the length of the subintervals beetwen t₀ and tₙ.
- `return_values::Bool`: If is true return t_values=[t₀,t₁...,tₙ],
   y_values=[y₀,y₁...,yₙ], else return yₙ.

## Example
```jldoctest
julia> using Qaylla

julia> F(t,y)=t-y
F (generic function with 1 method)

julia> euler(F,0,2,1,5)
0.9830400000000001

julia> euler(F,0,2,1,5,true)
(Any[0, 0.2, 0.4, 0.6000000000000001, 0.8, 1.0], Any[2, 1.6, 1.32, 1.1360000000000001, 1.0288000000000002, 0.9830400000000001])

julia> euler(F,0,2,1,0.2,true)
(Any[0, 0.2, 0.4, 0.6000000000000001, 0.8, 1.0], Any[2, 1.6, 1.32, 1.1360000000000001, 1.0288000000000002, 0.9830400000000001])
```
"""
function euler(F::Function,
               t₀::Number,
               y₀::Number,
               tₙ::Number,
               nh::Number,
               return_values::Bool=false)

    if nh>=1.0
        n=nh
        h=(tₙ-t₀)/n
    else
        n=(tₙ-t₀)/nh
        h=nh
    end

    n::Int64=n 
    yₙ=y₀

    if return_values
        t_values=[]
        y_values=[]

        append!(t_values,t₀)
        append!(y_values,y₀)

        tₙ=t₀

        for _ in 1:n
            yₙ=yₙ+h*F(tₙ,yₙ)
            tₙ=tₙ+h
            append!(y_values,yₙ)
            append!(t_values,tₙ)
        end

        return t_values,y_values
        
    end 


    for i in 1:n
        yₙ=yₙ+h*F(t₀+h*(i-1),yₙ)
    end

    return yₙ
end
