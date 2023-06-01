using WiltonInts84
using LinearAlgebra

# Equation numbers are referring to "Singularity extraction technique for integral equation methods
# with higher order basis functions on plane triangles and tetrahedra" by Järvenpää, Taskinen and Ylä-Oijala.
# Currently only for constant and linear basis functions, but should be possible to extend to higher
# order basis functions with equations given in the paper.
"""
function fundamentals(p,x)

    This function calculates all the needed properties of a triangle p, such as
    the normal, the lengths of the sides, the tangents of the sides etc.
"""
function fundamentals(p,x)
    T = eltype(p[1])
    e = Vector{Vector{T}}(undef,3)
    s = Vector{Vector{T}}(undef,3)
    m = Vector{Vector{T}}(undef,3)
    t0 = Vector{T}(undef,3)
    l = Vector{T}(undef,3)
    sm = Vector{T}(undef,3)
    sp = Vector{T}(undef,3)
    Rm = Vector{T}(undef,3)
    Rp = Vector{T}(undef,3)
    R0 = Vector{T}(undef,3)

    for i in 1:3
        e[i] = p[mod1(i+1,3)]-p[i]
    end
    n = normalize(cross(e[1],e[2]))
    for i in 1:3
        s[i] = normalize(e[mod1(i+1,3)])
        m[i] = cross(s[i],n)
    end
    û = s[3]
    v̂ = normalize(-s[2]+(s[2]⋅û)*û)
    ŵ = normalize(cross(û,v̂))

    u0 = (x-p[1])⋅û
    v0 = (x-p[1])⋅v̂
    w0 = (x-p[1])⋅ŵ

    for i in 1:3
        t0[i] = (p[mod1(i+1,3)]-x) ⋅ m[i]
        l[i] = norm(p[mod1(i+2,3)]-p[mod1(i+1,3)])
        sm[i] = (p[mod1(i+1,3)]-x)⋅s[i]
        sp[i] = sm[i] + l[i]
        Rp[i] = norm(x-p[mod1(i-1,3)])
        Rm[mod1(i+1,3)] = Rp[i]

        R0[i] = sqrt(t0[i]^2+w0^2)
    end
    u3 = û⋅(p[3]-p[1])
    v3 = v̂⋅(p[3]-p[1])
    l3 = l[3]
    return s, m, t0, l, sm, sp, Rm, Rp, R0, û, v̂, ŵ, u0,v0,w0,u3,v3,l3
end

""" 
function lineintegrals(m, l, sm, sp, Rm, Rp, R0, u3, v3, l3, N, T)

    Calculates integrals on the edges of a triangle with 
"""
function lineintegrals(m, l, sm, sp, Rm, Rp, R0,u3,v3,l3,N,T)
    # constant basis
    ∫∂ΔᵢR⁻¹ = zeros(3)
    for i in 1:3
        ∫∂ΔᵢR⁻¹[i] = log((sp[i]+Rp[i])/(sm[i]+Rm[i])) # eqn 31
    end
    ∫∂ΔᵢRⁿ = zeros(3,N÷2+3) # starting at R¹
    ∫∂ΔᵢRⁿ[:,1] = ∫∂ΔᵢR⁻¹
    for i in 1:3
        for j in 1:2:N+2
            ∫∂ΔᵢRⁿ[i,j÷2+2] = 1/(j+1)*(sp[i]*Rp[i]^j-sm[i]*Rm[i]^j+j*R0[i]^2*∫∂ΔᵢRⁿ[i,j÷2+1]) #eqn 33
        end
    end
    mᵢ∫∂ΔᵢRⁿ = Array{Vector{T}}(undef, size(∫∂ΔᵢRⁿ))
    for i in 1:3
        for j in 1:length(1:N÷2+3)
            mᵢ∫∂ΔᵢRⁿ[i,j] = m[i]*∫∂ΔᵢRⁿ[i,j]
        end
    end
    Σmᵢ∫∂ΔᵢRⁿ = sum(mᵢ∫∂ΔᵢRⁿ,dims = 1) # starting at R⁻¹

    # varying linearly on edge
    ∫∂ΔᵢRⁿs = zeros(3,N÷2+3) # starting at R⁻¹
    for i in 1:3
        for (z,j) in enumerate(-1:2:N+2)
            ∫∂ΔᵢRⁿs[i,z] = ((Rp[i]^(j+2)-Rm[i]^(j+2))/(j+2)) # eqn 36
        end
    end

    ∫∂ΔᵢRⁿu = ∫∂ΔᵢRⁿs.*[(u3-l3)/l[1],-u3/l[2],1] + ∫∂ΔᵢRⁿ.*[(l3*sp[1]-u3*sm[1])/l[1],u3*sp[2]/l[2],-sm[3]] # starting at R⁻¹, follows from eqn 36 and eqn 52
    ∫∂ΔᵢRⁿv = ∫∂ΔᵢRⁿs.*[v3/l[1],-v3/l[2],0] + ∫∂ΔᵢRⁿ.*[(-v3*sm[1])/l[1],v3*sp[2]/l[2],0] # starting at R⁻¹, follows from eqn 36 and eqn 52

    mᵢ∫∂ΔᵢRⁿu = Array{Vector{T}}(undef, size(∫∂ΔᵢRⁿu)) # starting at R⁻¹
    for i in 1:3
        for j in eachindex(∫∂ΔᵢRⁿu[1,:])
            mᵢ∫∂ΔᵢRⁿu[i,j] = m[i]*∫∂ΔᵢRⁿu[i,j]
        end
    end
    Σmᵢ∫∂ΔᵢRⁿu = sum(mᵢ∫∂ΔᵢRⁿu,dims = 1) # starting at R⁻¹

    mᵢ∫∂ΔᵢRⁿv = Array{Vector{T}}(undef, size(∫∂ΔᵢRⁿv)) # starting at R⁻¹
    for i in 1:3
        for j in eachindex(∫∂ΔᵢRⁿv[1,:])
            mᵢ∫∂ΔᵢRⁿv[i,j] = m[i]*∫∂ΔᵢRⁿv[i,j]
        end
    end
    Σmᵢ∫∂ΔᵢRⁿv = sum(mᵢ∫∂ΔᵢRⁿv,dims = 1) # starting at R⁻¹

    return ∫∂ΔᵢR⁻¹,∫∂ΔᵢRⁿ,Σmᵢ∫∂ΔᵢRⁿ,∫∂ΔᵢRⁿs,∫∂ΔᵢRⁿu,∫∂ΔᵢRⁿv,Σmᵢ∫∂ΔᵢRⁿu,Σmᵢ∫∂ΔᵢRⁿv
end


function higherorder(p1,p2,p3,x,N)
    T = eltype(p1)
    p = [p1,p2,p3]
    s, m, t0, l, sm, sp, Rm, Rp, R0, û, v̂, ŵ, u0,v0,w0,u3,v3,l3 = fundamentals(p,x)

    ∫∂ΔᵢR⁻¹,∫∂ΔᵢRⁿ,Σmᵢ∫∂ΔᵢRⁿ,∫∂ΔᵢRⁿs,∫∂ΔᵢRⁿu,∫∂ΔᵢRⁿv,Σmᵢ∫∂ΔᵢRⁿu,Σmᵢ∫∂ΔᵢRⁿv = lineintegrals(m, l, sm, sp, Rm, Rp, R0,u3,v3,l3,N,T)

    β = zeros(3)
    for i in 1:3
        β[i] = atan(t0[i]*sp[i]/(R0[i]^2+abs(w0)*Rp[i]))-atan(t0[i]*sm[i]/(R0[i]^2+abs(w0)*Rm[i])) # eqn 40
    end

    ∫ΔRⁿ = zeros(N÷2+3) # starting at R⁻³
    ∫ΔRⁿ[1] = w0==0.0 ? 0.0 : 1/abs(w0)*sum(β) # eqn 39
    for (j,i) in enumerate(-1:2:N) 
            ∫ΔRⁿ[j+1] = sum(t0./(i+2).*∫∂ΔᵢRⁿ[:,j]) + i * w0^2/(i+2)*∫ΔRⁿ[j] # eqn 44
    end

    ∫ΔRⁿu = zeros(N÷2+3) # starting at R⁻³
    for (i,j) in enumerate(-3:2:N+1)
            ∫ΔRⁿu[i] = û/(j+2)⋅Σmᵢ∫∂ΔᵢRⁿ[i] + u0*∫ΔRⁿ[i] #  eqn. 53
    end

    ∫ΔRⁿv = zeros(N÷2+3) # starting at R⁻³
    for (i,j) in enumerate(-3:2:N+1)
        ∫ΔRⁿv[i] = v̂/(j+2)⋅Σmᵢ∫∂ΔᵢRⁿ[i] + v0*∫ΔRⁿ[i] #  eqn. 54
    end

    linearbasis = [1 -1/l3 (u3/l3-1)/v3; 0 1/l3 -u3/l3/v3; 0 0 1/v3 ] # cdot [1 u v] basis as defined in e.g. Graglia 1993
    ∫NᵢRⁿ = Vector{Vector{T}}(undef,N÷2+3)
    for (i,j) in enumerate(-3:2:N+1)
        ∫NᵢRⁿ[i] = linearbasis * [∫ΔRⁿ[i], ∫ΔRⁿu[i], ∫ΔRⁿv[i]]
    end

    ∫Δ∇Rⁿu = Vector{Vector{T}}(undef, N÷2+2)
    for (i,j) in enumerate(-1:2:N+1)
        ∫Δ∇Rⁿu[i] = w0*ŵ*(û⋅Σmᵢ∫∂ΔᵢRⁿ[i]+u0*j*∫ΔRⁿ[i])-û*(û⋅Σmᵢ∫∂ΔᵢRⁿu[i])-v̂*(v̂⋅Σmᵢ∫∂ΔᵢRⁿu[i])+∫ΔRⁿ[i+1]*û # follows from eqn 56 and eqn 53/54
    end
    ∫Δ∇Rⁿv = Vector{Vector{T}}(undef, N÷2+2)
    for (i,j) in enumerate(-1:2:N+1)
        ∫Δ∇Rⁿv[i] = w0*ŵ*(v̂⋅Σmᵢ∫∂ΔᵢRⁿ[i]+v0*j*∫ΔRⁿ[i])-v̂*(v̂⋅Σmᵢ∫∂ΔᵢRⁿv[i])-û*(û⋅Σmᵢ∫∂ΔᵢRⁿv[i])+∫ΔRⁿ[i+1]*v̂ # follows from eqn 56 and eqn 53/54
    end
    ∫Δ∇Rⁿ = Vector{Vector{T}}(undef, N÷2+2)
    for (i,j) in enumerate(-1:2:N+1)
        ∫Δ∇Rⁿ[i] = -j*((∫ΔRⁿu[i]-u0*∫ΔRⁿ[i])*û+(∫ΔRⁿv[i]-v0*∫ΔRⁿ[i])*v̂ - w0*∫ΔRⁿ[i]*ŵ) # follows from eqn 56 and eqn 53/54
    end
    ∫Δ∇RⁿNᵢ = Vector{Vector{Vector{T}}}(undef,N÷2+2)
    for (i,j) in enumerate(-1:2:N+1)
        ∫Δ∇RⁿNᵢ[i] = linearbasis * [∫Δ∇Rⁿ[i], ∫Δ∇Rⁿu[i], ∫Δ∇Rⁿv[i]]
    end
 
    return ∫ΔRⁿ, ∫NᵢRⁿ,∫Δ∇Rⁿ, ∫Δ∇RⁿNᵢ#, ∫ΔRⁿu, ∫ΔRⁿv, ∫Δ∇Rⁿu, ∫Δ∇Rⁿv
end