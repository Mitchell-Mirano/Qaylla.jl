"""
    lagrange(x,y)

computes the lagrange interpolation polynomial

"""


function lagrange(x,y)::Function
    
    function p(value::Number)
        p_value=0
        n=length(x)

        for i in 1:n
            li=1
            for j in 1:n
                if i!=j
                    li*=(value-x[j])/(x[i]-x[j])
                end
            end
            p_value+=li*y[i]
        end
        
        return p_value
    end

    return p
end 
