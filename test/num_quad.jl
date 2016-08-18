using FastGaussQuadrature

function legendre(n,a,b)
    x,w = FastGaussQuadrature.gausslegendre(n)
    w *= (b-a)/2
    x = a+(x+1)*(b-a)/2
    return x, w
end

function dblquadints{N}(v1,v2,v3,x,T::Type{Val{N}})

    I = zeros(eltype(x),N+3)
    K = zeros(typeof(x),N+3)
    n = normalize(cross(v1-v3,v2-v3))
    ξ = x - dot(x-v1,n)*n

    # correctionn required for \int h/R^{-3} because
    # the definition of h is linked to triangle orientation
    σ = sign(dot(cross(v1-ξ,v2-ξ),n))
    J, L = dblquadints1(ξ,v1,v2,x,T)
    J *= σ; L *= σ;  J[1] *= σ
    I += J; K += L;

    σ = sign(dot(cross(v2-ξ,v3-ξ),n))
    J, L = dblquadints1(ξ,v2,v3,x,T)
    J *= σ;  L *= σ; J[1] *= σ
    I += J; K += L;

    σ = sign(dot(cross(v3-ξ,v1-ξ),n))
    J, L = dblquadints1(ξ,v3,v1,x,T)
    J *= σ; L *= σ; J[1] *= σ
    I += J; K += L;

    return I, K
end

function dblquadints1{N}(v1,v2,v3,x,::Type{Val{N}},ri=-1.0,ro=1.0e15)
    G = 100
    s, w = legendre(G, 0.0, 1.0)
    t1 = v1-v3
    t2 = v2-v3
    n = cross(t1,t2)
    I = zeros(eltype(x),N+3)
    K = zeros(typeof(x),N+3)
    a2 = norm(n)
    a2 < eps(eltype(x)) && return I, K
    n = normalize(n)
    d = dot(x-v1,n)
    ξ = x - d*n
    for g in 1:G
        u = s[g]
        for h in 1:G
            v = (1-u)*s[h]
            y = v3 + u*t1 + v*t2
            R = norm(x-y)
            (ri <= R <= ro) || continue
            j = a2 * w[g] * w[h] * (1-u)
            I[1] += j * d / R^3
            K[1] += j * (y-ξ) / R^3
            for n = -1 : N
                I[n+3] += j * R^n
                K[n+3] += j * (y-ξ) * R^n
            end
        end
    end
    return I, K
end
