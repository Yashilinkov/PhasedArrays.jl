"""
    abstract type ElementPattern end

Abstract supertype for representing the far-field radiation pattern of an antenna element.
Concrete subtypes should define how to evaluate the electric field components for given directions.
"""
abstract type ElementPattern end

"""
    IsotropicPattern(polarization::Char)

Represents an ideal isotropic radiator with a specified polarization.

# Arguments
- `polarization`: `'θ'` or `'ϕ'` indicating the direction of the non-zero electric field component.
"""
struct IsotropicPattern <: ElementPattern 
    polarization::Char
end

"""
    TabulatedPattern(theta, phi, Eθ, Eϕ)

Stores a tabulated element radiation pattern over a grid of `(θ, ϕ)` angles.

# Arguments
- `theta`: Vector of θ (theta) angles in degrees.
- `phi`: Vector of ϕ (phi) angles in degrees.
- `Eθ`: Matrix of complex electric field components in θ direction.
- `Eϕ`: Matrix of complex electric field components in ϕ direction.

# Notes
Assumes that `Eθ[i,j]` and `Eϕ[i,j]` correspond to `theta[i]` and `phi[j]`.
"""
struct TabulatedPattern <: ElementPattern
    theta::Vector{Float64}
    phi::Vector{Float64}
    Eθ::Matrix{ComplexF64}
    Eϕ::Matrix{ComplexF64}
end

"""
    element_gain(p::IsotropicPattern, θ, ϕ)

Returns the electric field components `(Eθ, Eϕ)` of an isotropic element
at a given direction `(θ, ϕ)`.

# Arguments
- `p`: An `IsotropicPattern` object.
- `θ`: Theta angle in degrees.
- `ϕ`: Phi angle in degrees.

# Returns
Tuple `(Eθ, Eϕ)` of complex values representing the element response.

# Notes
Returns unit magnitude for the specified polarization, and zero for the orthogonal component.
"""
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

"""
    element_gain(p::TabulatedPattern, θ, ϕ)

Looks up and returns the electric field components `(Eθ, Eϕ)` from the tabulated pattern.

# Arguments
- `p`: A `TabulatedPattern` object.
- `θ`: Theta angle in degrees (must match exactly one of the entries in `p.theta`).
- `ϕ`: Phi angle in degrees (must match exactly one of the entries in `p.phi`).

# Returns
Tuple `(Eθ, Eϕ)` of complex values at the specified direction.

# Warning
This function requires an exact match for `θ` and `ϕ` for now, interpolation is not yet implemented.
"""
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
