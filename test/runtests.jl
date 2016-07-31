using WiltonInts84
using Base.Test

# if Pkg.installed("FixedSizeArrays") == nothing
#     Pkg.add("FixedSizeArrays")
# end

using FixedSizeArrays
using FastGaussQuadrature

# write your own tests here
WI = WiltonInts84

v1 = Vec(0.0, 0.0, 0.0)
v2 = Vec(1.0, 0.0, 0.0)
v3 = Vec(0.0, 1.0, 0.0)
n = normalize(cross(v1-v3,v2-v3))
x = (v1 + v2 + v3)/3 + 20n

function legendre(n,a,b)
    x,w = FastGaussQuadrature.gausslegendre(n)
    w *= (b-a)/2
    x = a+(x+1)*(b-a)/2
    return x, w
end

function dblquadints{N}(v1,v2,v3,x,::Type{Val{N}})
    G = 20
    s, w = legendre(G, 0.0, 1.0)
    t1 = v1-v3
    t2 = v2-v3
    n = cross(t1,t2)
    a2 = norm(n)
    d = dot(x-v1,n)/a2
    I = zeros(eltype(x),N+3)
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

I = WI.wiltonints(v1,v2,v3,x,Val{7})
J = dblquadints(v1,v2,v3,x,Val{7})
maximum(abs(I-J)./abs(J)) < 1.0e-10

x = (1-1.5)*v1 + 1.5*v3 + 20*n
I = WI.wiltonints(v1,v2,v3,x,Val{7})
J = dblquadints(v1,v2,v3,x,Val{7})
