using WiltonInts84
using Base.Test

using FixedSizeArrays
using FastGaussQuadrature

# write your own tests here
WI = WiltonInts84

v1 = Vec(1.0, 0.0, 0.0)
v2 = Vec(0.0, 1.0, 0.0)
v3 = Vec(0.0, 0.0, 0.0)
n = normalize(cross(v1-v3,v2-v3))

function legendre(n,a,b)
    x,w = FastGaussQuadrature.gausslegendre(n)
    w *= (b-a)/2
    x = a+(x+1)*(b-a)/2
    return x, w
end

function dblquadints{N}(v1,v2,v3,x,T::Type{Val{N}})

    I = zeros(eltype(x),N+3)
    n = normalize(cross(v1-v3,v2-v3))
    ξ = x - dot(x-v1,n)*n

    σ = sign(dot(cross(v1-ξ,v2-ξ),n))
    K = σ * dblquadints1(ξ,v1,v2,x,T)
    # correctionn required for \int h/R^{-3} because
    # the definition of h is linked to triangle orientation
    K[1] *= σ; I += K

    σ = sign(dot(cross(v2-ξ,v3-ξ),n))
    K = σ * dblquadints1(ξ,v2,v3,x,T)
    K[1] *= σ; I += K

    σ = sign(dot(cross(v3-ξ,v1-ξ),n))
    K = σ * dblquadints1(ξ,v3,v1,x,T)
    K[1] *= σ; I += K

    return I
end

function dblquadints1{N}(v1,v2,v3,x,::Type{Val{N}})
    G = 20
    s, w = legendre(G, 0.0, 1.0)
    t1 = v1-v3
    t2 = v2-v3
    n = cross(t1,t2)
    I = zeros(eltype(x),N+3)
    a2 = norm(n)
    a2 < eps(eltype(x)) && return I
    d = dot(x-v1,n)/a2
    for g in 1:G
        u = s[g]
        for h in 1:G
            v = (1-u)*s[h]
            y = v3 + u*t1 + v*t2
            R = norm(x-y)
            j = a2 * w[g] * w[h] * (1-u)
            I[1] += j * d / R^3
            for n = -1 : N
                I[n+3] += j * R^n
            end
        end
    end
    return I
end

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

    (v1 + v2 + v3)/3               - 0n,  # h = 0, inside
    (1-1.5)v1 + 1.5v3              - 0n,  # h = 0, on extension [v1,v3]
    -0.5v1 + 0.5v2 + (1-0.5-0.5)v3 - 0n,  # h = 0, outside
]

τ = [
    1e-10,
    1e-10,
    1e-10,
    1e-10,
    1e-10,

    1e-10,
    1e-10,
    1e-10,
    1e-10,
    1e-10,

    1e-7,
    1e-10,
    1e-10,
]

function nearlyequal{T}(x::T,y::T,τ::T)
    X = abs(x)
    Y = abs(y)
    D = abs(x-y)

    x == y && return true

    ϵ = eps(zero(x))
    (x == 0 || y == 0 || D < ϵ) && return D < τ
    M = typemax(typeof(x))
    D / min(X+Y,M) < τ
end

for i in eachindex(X)
    x = X[i]
    I = WI.wiltonints(v1,v2,v3,x,Val{7})
    J = dblquadints(v1,v2,v3,x,Val{7})
    @test all(!isnan(I))
    @test all(!isinf(I))
    for j in eachindex(I)
        @test nearlyequal(I[j],J[j], τ[i])
    end
end


# # projection inside the triangle
# x = (v1 + v2 + v3)/3 + 20n
# I = WI.wiltonints(v1,v2,v3,x,Val{7})
# J = dblquadints(v1,v2,v3,x,Val{7})
# @test maximum(abs(I-J)./abs(J)) < 1.0e-10
#
# # projection outside of the triangle
# x = (v1 + v2 + v3)/3 - 20n
# I = WI.wiltonints(v1,v2,v3,x,Val{7})
# J = dblquadints(v1,v2,v3,x,Val{7})
# @test maximum(abs(I-J)./abs(J)) < 1.0e-10
#
# # projection on a triangle edge
# x = (1-1.5)*v1 + 1.5*v3 + 20*n
# I = WI.wiltonints(v1,v2,v3,x,Val{7})
# J = dblquadints(v1,v2,v3,x,Val{7})
# @test maximum(abs(I-J)./abs(J)) < 1.0e-10
#
# # projection outside of the triangle
# x = Vec(-0.5,0.5,20.0)
# I = WI.wiltonints(v1,v2,v3,x,Val{7})
# J = dblquadints(v1,v2,v3,x,Val{7})
# @test maximum(abs(I-J)./abs(J)) < 1.0e-10
#
# # projection on a vertex
# x = Vec(0.0,1.0,20.0)
# I = WI.wiltonints(v1,v2,v3,x,Val{7})
# J = dblquadints(v1,v2,v3,x,Val{7})
# @test maximum(abs(I-J)./abs(J)) < 1.0e-10
