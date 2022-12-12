module Yaqa

# interpolation methods
include("./interpolations/lagrange.jl")
include("./interpolations/newton.jl")

export lagrange,
       newton

# integral methods
include("./integrals/trapeze.jl")
include("./integrals/simpson.jl")
include("./integrals/newton-cotes.jl")
include("./integrals/romberg.jl")

export trapeze,
       simpson,
       newton_cotes,
       romberg

# integral2d methods
include("./2d-integrals/trapezed2.jl")
include("./2d-integrals/simpson2d.jl")

export trapeze2d,
       simpson2d

# differeential equations methods
include("./differential-equations/euler.jl")

export euler

end # module NumericalMethods
