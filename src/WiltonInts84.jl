module WiltonInts84

include("mutuple.jl")
include("contour.jl")


function segints!{N}(a, b, p, h, ::Type{Val{N}})

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


function arcints!{N}(α, m, p, h, ::Type{Val{N}})

    T = typeof(h)
    P = typeof(m)

    A = Vector{T}(N+3)
    B = Vector{T}(N+3)

    h2 = h^2
    p2 = p^2
    q2 = h2 + p2
    q = sqrt(q2)
    d = sqrt(h2)

    # n == -3
    A[1] = -α * (1/q - 1/d)
    B[1] = -p / q * m

    # n == -1
    A[2] = α * ((-A[1]/α + 1/d)*q2 - d)
    B[2] = p * q * m

    # n == 0
    A[3] = α * p2 / 2
    B[3] = p * q2 / 2 * m

    for i in 4:length(A)
        n = i-3
        A[i] = α * ((n*A[i-2]/α + d^n)*q2 - d^(n+2)) / (n+2)
        B[i] = p * q^(n+2) * m
    end

    A[1] *= h
    return A, B
end


function circleints!{N}(σ, p, h, ::Type{Val{N}})

    T = typeof(h)
    A = Vector{T}(N+3)

    h2 = h^2
    p2 = p^2
    q2 = h2 + p2
    q = sqrt(q2)
    d = sqrt(h2)
    α = σ * 2π

    # n == -3
    A[1] = -α * (1/q - 1/d)

    # n == -1
    A[2] = α * ((-A[1]/α + 1/d)*q2 - d)

    # n == 0
    A[3] = α * p2 / 2

    for i in 4:length(A)
        n = i-3
        A[i] = α * ((n*A[i-2]/α + d^n)*q2 - d^(n+2)) / (n+2)
    end

    A[1] *= h
    return A
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

        P,Q = segints!(sa, sb, p, h, UB)
        for i in eachindex(I) I[i] += P[i]   end
        for i in eachindex(K) K[i] += Q[i]*m end
    end

    return I, K
end



"""
    angle(p,q)

Returns the positive angle in [0,π] between p and q
"""
function angle(p,q)
    cs = dot(p,q) / norm(p) / norm(q)
    cs = clamp(cs, -1, +1)
    return acos(cs)
end


function wiltonints2{N}(ctr, x, UB::Type{Val{N}})

    n = ctr.normal
    h = ctr.height

    ξ = x - h*n

    I = zeros(eltype(x), N+3)
    K = zeros(typeof(x), N+3)

    # segments contributions
    for i in eachindex(ctr.segments)
        a = ctr.segments[i][1]
        b = ctr.segments[i][2]
        t = b - a
        t /= norm(t)
        m = cross(t, n)
        p = dot(a-ξ,m)
        sa = dot(a-ξ,t)
        sb = dot(b-ξ,t)
        P,Q = segints!(sa, sb, p, h, UB)
        for j in eachindex(I) I[j] += P[j]   end
        for j in eachindex(K) K[j] += Q[j]*m end
    end

    # arc contributions
    for i in eachindex(ctr.arcs)
        a = ctr.arcs[i][1]
        b = ctr.arcs[i][2]
        σ = crt.arcs[i][3]
        p = σ > 0 ? ctr.outer_radius : ctr.inner_radius
        u1 = (a - ξ) / p
        u2 = σ * (n × u1)
        ξb = b - ξ
        α = dot(ξb,u2) >= 0 ? σ*angle(ξb,u1) : σ*angle(ξb,-u1) + π
        m = σ*(sin(α)*u1 + (1-cos(α))*u2)
        P, Q = arcints!(α, m, p, h, UB)
        I .+= P
        K .+= Q
    end

    # circle contributions
    for i in eachindex(ctr.circles)
        σ = ctr.circles[i]
        p = σ > 0 ? ctr.outer_radius : ctr.inner_radius
        P = circleints!(σ, p, h, UB)
        I .+= P
    end

    return I, K
end



end # module
