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