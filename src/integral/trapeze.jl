function trapeze(f::Function,a::Number,b::Number,n::Int64=1)::Float64
    
    h=(b-a)/n
    ∑=0
    for i in 1:n-1 
        ∑+=f(a+h*i)
    end
    return (h/2)*(f(a)+2*∑+f(b))
    
end