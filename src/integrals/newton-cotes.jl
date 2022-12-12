"""
    newton_cotes(f::Function,a::Number,b::Number,n::Int64,closed::Bool=true)::Float64

Computes the integral of f applying the Newton-Cotes methods.
If closed=true `1≤n≤4`, else `0≤n≤3`.

## Example
```jldoctest
julia> using Qaylla

julia> f(x)=sin(x)
f (generic function with 1 method)

julia> a=0
0

julia> b=pi/2
1.5707963267948966

julia> newton_cotes(f,a,b,1)
0.7853981633974483

julia> newton_cotes(f,a,b,2)
1.0022798774922104

julia> newton_cotes(f,a,b,3)
1.001004923314279

julia> newton_cotes(f,a,b,4)
0.9999915654729927

julia> newton_cotes(f,a,b,0,false)
1.1107207345395915

julia> newton_cotes(f,a,b,1,false)
1.0728738432865557

julia> newton_cotes(f,a,b,2,false)
0.9979892924561773

julia> newton_cotes(f,a,b,3,false)
0.9986082958707345
```
"""
function newton_cotes(f::Function,a::Number,b::Number,n::Int64,closed::Bool=true)::Float64
    
    if closed==true
        h=(b-a)/n
        if n==1
            return (h/2)*(f(a)+ f(b))
        end
        if n==2
            return (h/3)*(f(a)+4*f(a+h)+f(b))
        end
        if n==3
            return (3*h/8)*(f(a)+3*f(a+h)+3*f(a+2*h)+f(b))
        end
        if n==4
            return (2*h/45)*(7*f(a)+32*f(a+h)+12*f(a+2*h)+32*f(a+3*h)+7*f(b))
        end
    end
    if closed==false
        h=(b-a)/(n+2)
        if n==0
            return (2*h)*f(a+h)
        end
        if n==1
            return (3*h/2)*(f(a+h)+f(a+2*h))
        end
        if n==2
            return (4*h/3)*(2*f(a+h)-f(a+2*h)+2*f(a+3*h))
        end
        if n==3
            return (5*h/24)*(11*f(a+h)+f(a+2*h)+f(a+3*h)+11*f(a+4*h))
        end
    end

end