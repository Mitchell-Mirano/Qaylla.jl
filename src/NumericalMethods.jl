module NumericalMethods

    # interpolation methods
    include("./interpolate/lagrange.jl")
    include("interpolate/newton.jl")

    export lagrange
           newton

    # integral methods
    include("./integral/trapeze.jl")
    include("./integral/simpson.jl")
    include("./integral/newton-cotes.jl")
    include("./integral/romberg.jl")
    
    export trapeze
           simpson
           newtonCotes
           romberg
end # module NumericalMethods
