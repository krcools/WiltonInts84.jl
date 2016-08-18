using WiltonInts84
WI = WiltonInts84

using Base.Test
using FixedSizeArrays

include("num_quad.jl")


v1 = Vec(1.0, 0.0, 0.0)
v2 = Vec(0.0, 1.0, 0.0)
v3 = Vec(0.0, 0.0, 0.0)
n = normalize(cross(v1-v3,v2-v3))



X = [
    (v1 + v2 + v3)/3               + 20n, # h > 0, inside
    (1-1.5)v1 + 1.5v3              + 20n, # h > 0, on extension [v1,v3]
    0.5v1 + 0.5v3                  + 20n, # h > 0, on the interior of [v1,v3]
    -0.5v1 + 0.5v2 + (1-0.5-0.5)v3 + 20n, # h > 0, outside
    v2                             + 20n, # h > 0, on top of v2

    (v1 + v2 + v3)/3               - 20n, # h < 0, inside
    (1-1.5)v1 + 1.5v3              - 20n, # h < 0, on extension [v1,v3]
    0.5v1 + 0.5v3                  - 20n, # h < 0, on the interior of [v1,v3]
    -0.5v1 + 0.5v2 + (1-0.5-0.5)v3 - 20n, # h < 0, outside
    v2                             - 20n, # h < 0, on top of v2
]

Y = [
    (1-1.5)v1 + 1.5v3              - 0n,  # h = 0, on extension [v1,v3]
    -0.5v1 + 0.5v2 + (1-0.5-0.5)v3 - 0n,  # h = 0, outside
]

Z = [
    (v1 + v2 + v3)/3               - 0n,  # h = 0, inside
]


function nearlyequal{T,U}(x::T,y::T,τ::U)

    X = norm(x)
    Y = norm(y)
    D = norm(x-y)

    x == y && return true

    ϵ = eps(one(τ))
    (X <= τ && Y <= τ) && return true
    2 * D / (X+Y) < τ
end

for (i,x) in enumerate(X)
    ctr = WI.contour(v1,v2,v3,x,-1.0,100.0)
    I, K = WI.wiltonints2(ctr,x,Val{7})
    #I, K = WI.wiltonints(v1,v2,v3,x,Val{7})
    J, L = dblquadints(v1,v2,v3,x,Val{7})
    @test all(!isnan(I))
    @test all(!isinf(I))
    for j in eachindex(I)
        @test nearlyequal(I[j],J[j], 1.0e-6)
        @test nearlyequal(K[j],L[j],1.0e-6)
    end
end

# Tests where a singularity cancelation
# approach is unfortunately not possible
for (i,x) in enumerate(Y)
    ctr = WI.contour(v1,v2,v3,x,-1.0,100.0)
    I, K = WI.wiltonints2(ctr,x,Val{7})
    #I, K = WI.wiltonints(v1,v2,v3,x,Val{7})
    J, L = dblquadints1(v1,v2,v3,x,Val{7})
    @test all(!isnan(I))
    @test all(!isinf(I))
    for j in eachindex(I)
        @test nearlyequal(I[j],J[j], 1.0e-6)
        @test nearlyequal(K[j],L[j],1.0e-6)
    end
end

# Test that only make sense for the scalar ints
for x in Z
    ctr = WI.contour(v1,v2,v3,x,-1.0,100.0)
    I, K = WI.wiltonints2(ctr,x,Val{7})
    #I, K = WI.wiltonints(v1,v2,v3,x,Val{7})
    J, L = dblquadints(v1,v2,v3,x,Val{7})
    @test all(!isnan(I))
    @test all(!isinf(I))
    for j in eachindex(I)
        @test nearlyequal(I[j],J[j], 1.0e-6)
    end
end
