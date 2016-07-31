module WiltonInts84

include("mutuple.jl")


"""
    oddints!(I, a, b, p, h Val{N})

computes the integrals

```math
I_{n} = \frac{p}{n+2} \int_a^b \frac{(s^2+p^2+h^2)^{(n+2)/2}}{s^2+p^2} ds
```

for ``n`` in `-3:2:N`. `N` is required to be odd.

**note**: The first element of the returned array is in fact ``h I_{-3}``. This
is done to avoid singularities when ``h=0``.
"""
function oddints!{N}(I, a, b, p, h, ::Type{Val{N}})

    @assert N ≥ -3
    @assert length(I) == N+3

    q2 = p^2+h^2; q = √q2
    K = zeros(I)
    K[1] += sign(h) * atan( (p*b) / (q2 + abs(h)*√(b^2+q2)) )
    K[1] -= sign(h) * atan( (p*a) / (q2 + abs(h)*√(a^2+q2)) )
    J = log(b/q + √(1+(b/q)^2)) - log(a/q + √(1+(a/q)^2))
    for i in 2 : 2 : length(I)
        n = i-3
        Kp = (i == 2) ? K[1]/h : K[i-2]
        # Compute I_{n} using J_{n} and I_{n+1}
        K[i] += p/(n+2)*(q)^(n+1)*J + n/(n+2)*h^2*Kp
        # Now compute J_{n+2}
        J = (b/q*√(1+(b/q)^2)^(n+2) - a/q*√(1+(a/q)^2)^(n+2) + (n+2)*J) / (n+3)
    end

    for i in eachindex(K) I[i] += K[i] end

end


function evenints!{N}(I, a, b, p, h, UB::Type{Val{N}})

    @assert N ≥ -3
    @assert length(I) == N+3
    M = div(N,2)

    pa = zeros(typeof(a), N+2); pa[1] = 1;
    pb = zeros(typeof(a), N+2); pb[1] = 1;
    pp = zeros(typeof(a), M+1); pp[1] = 1;
    hp = zeros(typeof(a), M+1); hp[1] = 1;

    # compute the required powers of a, b, and p
    for i in 2:length(pa) pa[i] = a * pa[i-1] end
    for i in 2:length(pb) pb[i] = b * pb[i-1] end
    for i in 2:length(pp) pp[i] = p*p * pp[i-1] end
    for i in 2:length(pp) hp[i] = h*h * hp[i-1] end

    # get the binomial coeffs
    bn = binomial(N+1)

    # compute the degree independent coefficients (see notebook for details)
    Ca = ckps(pp,pa,bn,UB)
    Cb = ckps(pp,pb,bn,UB)

    # now put the actual integrals together
    for i in 3 : 2 : length(I)
        n = i-3
        @assert iseven(n)
        m = div(n,2)
        α = p/(n+2)
        for k = 1 : m+1
            I[i] += α * bn[k+1,m+2] * hp[m-k+2] * (Cb[k] - Ca[k])
        end
    end


    return I
end


function ckps{N}(p,s,bn,::Type{Val{N}})


    m = div(N,2)
    C = zeros(eltype(p), m+1)
    for k in 1 : m+1
        for l in 0 : k-1
            F1 = bn[l+1,k]
            F2 = s[2l+2]
            F3 = p[k-l]
            C[k] += bn[l+1,k]*s[2l+2]*p[k-l]/(2l+1)
        end
    end

    return C
end


"""
    binomial(N)

Compute the table of binomial coefficients up to (N+1,.)
"""
function binomial(N)
    bn = zeros(Int, N+2, N+2)
    for n in 0:N+1 (j = n+1; bn[1,j] = 1) end
    for n in 1:N+1
        j = n+1
        for k in 1:n
            i = k+1
            bn[i,j] = bn[i-1,j-1] + bn[i,j-1]
        end
    end
    return bn
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

        oddints!(I, sa, sb, p, h, UB)
        evenints!(I, sa, sb, p, h, UB)
    end

    return I
end


end # module
