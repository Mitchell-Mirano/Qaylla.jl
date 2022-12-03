function euler(f::Function,
               t0::Number,
               y0::Number,
               tn::Number,
               nh::Number,
               return_values::Bool=false)


    if nh>=1.0
        n=nh
        h=(tn-t0)/n
    else
        n=(tn-t0)/nh
        h=nh
    end

    n::Int64=n 
    yn=y0 

    if return_values
        t_values=[]
        y_values=[]

        append!(t_values,t0)
        append!(y_values,y0)

        tn=t0

        for _ in 1:n
            yn=yn+h*f(tn,yn)
            tn=tn+h
            append!(y_values,yn)
            append!(t_values,tn)
        end

        return t_values,y_values
        
    end 


    for i in 1:n
        yn=yn+h*f(t0+h*(i-1),yn)
    end

    return yn
end

