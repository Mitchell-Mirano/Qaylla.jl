function romberg(f::Function,a::Number,b::Number,n::Int64=4)::Float64
    
    h=(b-a)
    R=zeros(n,n)
    R[1,1]=(h/2)*(f(a)+f(b))
    for k in 2:n
        ∑=0
        for i in 1:2^(k-2)
            ∑+=f(a+(i-0.5)*h)
        end
        R[k,1]=(0.5)*(R[k-1,1]+h*∑)
        h=h/2
    end

    for i in 2:n
        for j in 2:n
            if j<=i
                R[i,j]=(4^(j-1)*R[i,j-1]-R[i-1,j-1])/(4^(j-1)-1)
            end
        end
    end
    
    return R[n,n]
    
end

