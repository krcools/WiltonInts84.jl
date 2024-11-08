using LinearAlgebra

include("num_quad.jl")

include("test_triangleints.jl")
include("test_segintsg.jl")
include("test_contour0.jl")
include("test_contour1.jl")
include("test_contour2.jl")
include("scan1.jl")
include("issue1.jl")
include("test_higherorderints.jl")

using TestItemRunner
@run_package_tests