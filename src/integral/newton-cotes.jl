function newtonCotes(f::Function,a::Number,b::Number,n::Int64,closed::Bool=true)::Float64
    
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