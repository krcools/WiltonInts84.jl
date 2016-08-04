module WiltonInts84

include("mutuple.jl")


function allints!{N}(a, b, p, h, ::Type{Val{N}})

    @assert N ≥ 0
    T = typeof(a)

    h2, d = h^2, abs(h)
    q2 = p^2+h2
    ra, rb = sqrt(a^2+q2), sqrt(b^2+q2)

    I = zeros(T, N+3)
    K = zeros(T, N+3)
    J = zeros(T, 2)

    # n = -3
    if p == 0
        I[1] = 0
    else
        I[1] = sign(h)*(atan((p*b)/(q2+d*rb)) - atan((p*a)/(q2 + d*ra)))
    end # I_{-3}
    if q2 == zero(typeof(q2))
        J[1] = (b > 0) ? log(b/a) : log(a/b)
    else
        J[1] = log(b + rb) - log(a + ra)
    end # J_{-1}
    K[1] = -J[1] # K_{-3}

    # n = -1
    I[2] = p*J[1] - h*I[1]
    J[1] = (b*rb - a*ra + q2*J[1])/2 # J_1
    K[2] = J[1]

    # n = 0
    I[3] = (b*p - a*p)/2
    J[2] = ((b*(b^2+q2)+2*q2*b) - (a*(a^2+q2)+2*q2*a))/3 # J_2
    K[3] = J[2]/2

    # n >= 1
    for i in 4 : length(I)
        n = i - 3; j = mod1(n,2)
        I[i] = p/(n+2)*J[j] + n/(n+2)*h2*I[i-2]
        J[j] = (b*rb^(n+2) - a*ra^(n+2) + (n+2)*q2*J[j])/(n+3)
        K[i] = J[j]/(n+2)
    end

    return I, K
end


function wiltonints{N}(v1,v2,v3,x,UB::Type{Val{N}})

    n = cross(v1-v3, v2-v3)
    n /= norm(n)

    h = dot(x-v1,n)
    ξ = x - h*n
    inside = true

    I = zeros(eltype(x), N+3)
    K = zeros(typeof(x), N+3)

    # compute the contributions from the edges of the triangle
    Va = (v1,v2,v3)
    Vb = (v2,v3,v1)
    for i in 1:3
        va = Va[i]
        vb = Vb[i]

        t = vb - va
        t /= norm(t)
        m = cross(t, n)
        p = dot(va-ξ,m)
        if p < 0
            inside = false
        end
        sa = dot(va-ξ,t)
        sb = dot(vb-ξ,t)

        P,Q = allints!(sa, sb, p, h, UB)
        for i in eachindex(I) I[i] += P[i]   end
        for i in eachindex(K) K[i] += Q[i]*m end
    end

    return I, K
end


end # module
