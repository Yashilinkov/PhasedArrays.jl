abstract type ElementPattern end

struct IsotropicPattern <: ElementPattern 
    polarization::Char
end

struct TabulatedPattern <: ElementPattern
    theta::Vector{Float64}
    phi::Vector{Float64}
    Eθ::Matrix{ComplexF64}
    Eϕ::Matrix{ComplexF64}
end

function element_gain(p::IsotropicPattern, θ, ϕ)
    if p.polarization == 'θ'
        Eθ = 1.0 + 0.0im
        Eϕ = 0.0 + 0.0im
    elseif p.polarization == 'ϕ'
        Eθ = 0.0 + 0.0im
        Eϕ = 1.0 + 0.0im
    end
    return Eθ,Eϕ 
end

function element_gain(p::TabulatedPattern, θ, ϕ)
    i = findfirst(p.theta .== θ)
    j = findfirst(p.phi .== ϕ)
    
    return p.Eθ[i,j], p.Eϕ[i,j]
end

function compute_dipole_pattern(axis::Char, θ, ϕ)
 
    if axis == 'z'
        Eθ = sind(θ) + 0.0im
        Eϕ = 0.0 + 0.0im
    elseif axis == 'y'
        Eθ =  cosd(θ) * sind(ϕ) + 0.0im
        Eϕ =  cosd(ϕ) + 0.0im
    elseif axis == 'x'
        Eϕ =  cosd(θ) * cosd(ϕ) + 0.0im
        Eθ =  sind(ϕ) + 0.0im
    end


    return Eθ,Eϕ
end
