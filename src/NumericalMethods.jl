module NumericalMethods

    include("./integral/trapeze.jl")
    include("./integral/simpson.jl")
    include("./integral/newton-cotes.jl")
    include("./integral/romberg.jl")
    
    export trapeze
           simpson
           newtonCotes
           romberg
end # module NumericalMethods
