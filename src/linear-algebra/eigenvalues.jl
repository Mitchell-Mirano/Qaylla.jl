using LinearAlgebra

function max_eigenvals(A,x,n)
    x=transpose(x)
    for _ in 1:n
        y=A*x
        x=y/norm(y)
    end
    λ=transpose(transpose(x)*A*x)[1]
    return λ,transpose(x)
end

function min_eigenvals(A,x,n)
    x=transpose(x)
    for _ in 1:n
        y=inv(A)*x
        x=y/norm(y)
    end
    λ=transpose(transpose(x)*A*x)[1]
    return λ,transpose(x)
end

A=[2 1 2;
   1 4 1; 
   1 2 2 ]

x=[1 0 1]

λ,x=min_eigenvals(A,x,10)
println("λ=$λ")
println("x=$x")