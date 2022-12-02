function simpson(f::Function,a::Number,b::Number,n::Int64=2)::Float64
    
    if n%2!=0
        n=n+1
    end
    h=(b-a)/n
    m=Integer(n/2)

    ∑1=0
    for i in 1:m
        ∑1+=f(a+(2*i-1)*h)
    end

    ∑2=0
    for i in 1:m-1
        ∑2+=f(a+(2*i)*h)
    end

    return (h/3)*(f(a)+f(b)+4*∑1+2*∑2)
    
end