using WiltonInts84
using Test

using StaticArrays
using LinearAlgebra

for U in [Float32, Float64]
    x1 = SVector{3,U}(1.0, 0.0, 0.0)
    x2 = SVector{3,U}(0.0, 0.0, 0.0)
    x3 = SVector{3,U}(0.5, sqrt(3)/2, 0.0)
    x = SVector{3,U}(-0.1, -sqrt(1e-8), 0.0)


    ctr=contour(x1, x2, x3, x)
    n = ctr.normal
    h = ctr.height
    ξ = x - h*n
    UB=Val{1}
    for j in 1:3
        a = ctr.segments[j][1]
        b = ctr.segments[j][2]
        t = b - a
        t /= norm(t)
        m = cross(t, n)
        p = dot(a-ξ, m)
        sa = dot(a-ξ, t)
        sb = dot(b-ξ, t)
        P, Q = WiltonInts84.segintsg(sa, sb, p, h, m, UB)
        for i in P
            @test !isnan(i)
            @test !isinf(i)
        end
    end
end


