module WiltonInts84

include("mutuple.jl")


function allints!{N}(I, a, b, p, h, ::Type{Val{N}})

    @assert N ≥ 3
    @assert length(I) == N+3

    h2 = h^2
    q2 = p^2+h2
    q = sqrt(q2)

    K = zeros(I)
    J = zeros(eltype(I),2)

    # n = -3
    if p == 0
        K[1] = 0
    else
        K[1] += sign(h) * atan( (p*b) / (q2 + abs(h)*√(b^2+q2)) )
        K[1] -= sign(h) * atan( (p*a) / (q2 + abs(h)*√(a^2+q2)) )
    end

    if q2 == zero(typeof(q2))
        J[1] = (b > 0) ? log(b/a) : log(a/b)
    else
        J[1] = log(b + sqrt(b^2 + q2)) - log(a + sqrt(a^2 + q2))
    end

    # n = -1
    K[2] = p*J[1] - h*K[1]
    J[1] = (b*√(b^2+q2) - a*√(a^2+q2) + q2*J[1])/2

    # n = 0
    K[3] += 0.5*b*p
    K[3] -= 0.5*a*p
    J[2] = ((b*(b^2+q2)+2*q2*b) - (a*(a^2+q2)+2*q2*a))/3

    # n >= 1
    for i in 4 : length(K)
        n = i - 3; j = mod1(n,2)
        K[i] = p/(n+2)*J[j] + n/(n+2)*h2*K[i-2]
        J[j] = (b*√(b^2+q2)^(n+2) - a*√(a^2+q2)^(n+2) + (n+2)*q2*J[j])/(n+3)
    end

    for i in eachindex(K) I[i] += K[i] end
end


function wiltonints{N}(v1,v2,v3,x,UB::Type{Val{N}})

    n = cross(v1-v3, v2-v3)
    n /= norm(n)

    h = dot(x-v1,n)
    ξ = x - h*n
    inside = true

    I = zeros(eltype(x), N+3)

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

        allints!(I, sa, sb, p, h, UB)
    end

    return I
end


end # module
