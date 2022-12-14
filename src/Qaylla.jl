module Qaylla

# interpolation methods
include("./interpolations/lagrange.jl")
include("./interpolations/newton.jl")
include("./interpolations/cubic-splines.jl")

export lagrange,
       newton,
       cubic_splines

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
include("./differential-equations/modified-euler.jl")
include("differential-equations/runge-kutta2.jl")
include("differential-equations/heun.jl")
include("differential-equations/runge-kutta4.jl")

export euler, 
       modified_euler,
       runge_kutta_order_2,
       heun,
       runge_kutta_order_4

end # module NumericalMethods
