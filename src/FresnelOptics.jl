module FresnelOptics

const VecI  = AbstractVector # Input  Vector
const VecO  = AbstractVector # Output Vector
const VecB  = AbstractVector # Buffer Vector
const VecIO = AbstractVector # In/Out Vector
const MatI  = AbstractMatrix # Input  Matrix
const MatO  = AbstractMatrix # Output Matrix
const MatB  = AbstractMatrix # Buffer Matrix
const MatIO = AbstractMatrix # In/Out Matrix

const 𝚷 = 2.0 * π

using LinearAlgebra.BLAS: gemv!

include("./utils/la.jl")
include("./fresnel.jl")
include("./lorentz.jl")
include("./dielectric.jl")
include("./layers.jl")

end # module
