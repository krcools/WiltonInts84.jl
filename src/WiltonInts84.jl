module WiltonInts84

export contour, wiltonints

include("contour.jl")

function segints!{N}(P, Q, a, b, p, h, m, ::Type{Val{N}})

    @assert N ≥ 0
    T = typeof(a)
    z = zero(T)

    h2, d = h^2, abs(h)
    q2 = p^2+h2
    ra, rb = sqrt(a^2+q2), sqrt(b^2+q2)

    I = zeros(T, N+3)
    K = zeros(T, N+3)
    J = (z,z)

    # n = -3
    sgn = norm(h) < eps(T)*1e3 ? zero(T) : sign(h)
    I[1] = p == 0 ? 0 : sgn*(atan((p*b)/(q2+d*rb)) - atan((p*a)/(q2 + d*ra)))
    j = q2 == 0 ? (b > 0 ? log(b/a) : log(a/b)) : log(b + rb) - log(a + ra)
    J = (j,z)
    K[1] = -j

    j = z
    J = (j,J[1])

    # n = -1
    I[2] = p*J[2] - h*I[1]
    j = (b*rb - a*ra + q2*J[2])/2
    J = (j,J[1])
    K[2] = j

    # n = 0
    I[3] = (b*p - a*p)/2
    j = ((b*(b^2+q2)+2*q2*b) - (a*(a^2+q2)+2*q2*a))/3
    J = (j,J[1])
    K[3] = j/2

    # n >= 1
    for i in 4 : length(I)
        n = i - 3 #; j = mod1(n,2)
        I[i] = p/(n+2)*J[2] + n/(n+2)*h2*I[i-2]
        j = (b*rb^(n+2) - a*ra^(n+2) + (n+2)*q2*J[2])/(n+3)
        J = (j,J[1])
        K[i] = j/(n+2)
    end

    for k in eachindex(P) P[k] += I[k]   end
    for k in eachindex(Q) Q[k] += K[k]*m end
end


@generated function segints!{N}(P, Q, a, b, p, h, m, ::Type{Val{N}})

  @assert N >= 0
  N3 = N +3

  xp = quote

    h2, d = h^2, abs(h)
    q2 = p^2+h2
    ra, rb = sqrt(a^2+q2), sqrt(b^2+q2)

    # n = -3
    sgn = norm(h) < eps(T)*1e3 ? zero(T) : sign(h)
    I1 = p == 0 ? 0 : sgn*(atan((p*b)/(q2+d*rb)) - atan((p*a)/(q2 + d*ra)))
    j = q2 == 0 ? (b > 0 ? log(b/a) : log(a/b)) : log(b + rb) - log(a + ra)
    J = (j,z)
    K1 = -j

    j = z
    J = (j,J1)

    # n = -1
    I2 = p*J[2] - h*I1
    j = (b*rb - a*ra + q2*J[2])/2
    J = (j,J[1])
    K2 = j

    # n = 0
    I3 = (b*p - a*p)/2
    j = ((b*(b^2+q2)+2*q2*b) - (a*(a^2+q2)+2*q2*a))/3
    J = (j,J[1])
    K3 = j/2

  end

  for i in 4 : N3
    Ip = symbol(:I,i-2)
    In = symbol(:I,i)
    Kn = symbol(:K,i)
    it = quote
      n = i - 3
      $In = p/(n+2)*J[2] + n/(n+2)*h2*$Ip
      j = (b*rb^(n+2) - a*ra^(n+2) + (n+2)*q2*J[2])/(n+3)
      J = (j,J[1])
      $Kn = j/(n+2)
    end
    append!(xp.args, it.args)
  end

  xpi = Expr(:tuple)
  xpk = Expr(:tuple)
  for i in 1 : N3
    push!(xpi.args, symbol(:I,i))
    push!(xpk.args, symbol(:K,i))
  end

  push!(xp.args, :($xpi, $xpk))

  return xp
end


function arcints!{N}(I, K, α, m, p, h, ::Type{Val{N}})

    T = typeof(h)
    P = typeof(m)

    A = Vector{T}(N+3)
    B = Vector{P}(N+3)

    h2 = h^2
    p2 = p^2
    q2 = h2 + p2
    q = sqrt(q2)
    d = sqrt(h2)

    # n == -3
    sgn = norm(h) < eps(T)*1e3 ? zero(T) : sign(h)
    A[1] = -α * (h/q - sgn)
    B[1] = -p / q * m

    # n == -1
    A[2] = α * (q - d)
    B[2] = p * q * m

    # n == 0
    A[3] = α * p2 / 2
    B[3] = p * q2 / 2 * m

    for i in 4:length(A)
        n = i-3
        A[i] = α * ((n*A[i-2]/α + d^n)*q2 - d^(n+2)) / (n+2)
        B[i] = p * q^(n+2) * m / (n+2)
    end

    for j in eachindex(I) I[j] += A[j]   end
    for j in eachindex(K) K[j] += B[j]   end
end


function circleints!{N}(I, K, σ, p, h, ::Type{Val{N}})

    T = typeof(h)
    A = Vector{T}(N+3)

    d = norm(h)
    h2 = h^2
    p2 = p^2
    q2 = h2 + p2
    q = sqrt(q2)
    α = σ * 2π

    # n == -3
    sgn = norm(h) < eps(T)*1e3 ? zero(T) : sign(h)
    A[1] = -α * (h/q - sgn)

    # n == -1
    A[2] = α * (q - d)

    # n == 0
    A[3] = α * p2 / 2

    for i in 4:length(A)
        n = i-3
        A[i] = α * ((n*A[i-2]/α + d^n)*q2 - d^(n+2)) / (n+2)
    end

    for j in eachindex(I) I[j] += A[j]   end
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


function wiltonints{N}(ctr, x, UB::Type{Val{N}})

    n = ctr.normal
    h = ctr.height

    ξ = x - h*n

    I = zeros(eltype(x), N+3)
    K = zeros(typeof(x), N+3)
    G = zeros(typeof(x), N+3)

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
        segints!(I, K, sa, sb, p, h, m, UB)
    end

    # arc contributions
    for i in eachindex(ctr.arcs)
        a = ctr.arcs[i][1]
        b = ctr.arcs[i][2]
        σ = ctr.arcs[i][3]
        p = σ > 0 ? ctr.plane_outer_radius : ctr.plane_inner_radius
        u1 = (a - ξ) / p
        u2 = σ * (n × u1)
        ξb = b - ξ
        α = dot(ξb,u2) >= 0 ? σ*angle(ξb,u1) : σ*(angle(ξb,-u1) + π)
        m = (sin(α)*u1 + σ*(1-cos(α))*u2)
        arcints!(I, K, α, m, p, h, UB)
    end

    # circle contributions
    for i in eachindex(ctr.circles)
        σ = ctr.circles[i]
        p = σ > 0 ? ctr.plane_outer_radius : ctr.plane_inner_radius
        circleints!(I, K, σ, p, h, UB)
    end

    G[1] = K[1] - I[1]
    for i in 2:length(G)
        G[i] = K[i] - h*I[i]
    end

    return I, K, G
end

"""
    wiltonints(p1,p2,p3,x,[r,R],Val{N})

Compute potential integrals over a triangle (intersected with a spherical)
mesh. Powers of the distance up to degree `N` are computed.
"""
function wiltonints{N}(p1,p2,p3,x,r,R,VN::Type{Val{N}})
    ctr = contour(p1,p2,p3,x,r,R)
    wiltonints(ctr,x,VN)
end

function wiltonints{N}(p1,p2,p3,x,VN::Type{Val{N}})
    ctr = contour(p1,p2,p3,x)
    wiltonints(ctr,x,VN)
end




end # module
