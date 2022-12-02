module NumericalMethods

       # interpolation methods
       include("./interpolate/lagrange.jl")
       include("./interpolate/newton.jl")

       export lagrange,
              newton

       # integral methods
       include("./integral/trapeze.jl")
       include("./integral/simpson.jl")
       include("./integral/newton-cotes.jl")
       include("./integral/romberg.jl")
       
       export trapeze,
              simpson,
              newtonCotes,
              romberg

       # integral2d methods
       include("./integral2d/trapezed2.jl")
       include("./integral2d/simpson2d.jl")

       export trapeze2d,
              simpson2d

end # module NumericalMethods
