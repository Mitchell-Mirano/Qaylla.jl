
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