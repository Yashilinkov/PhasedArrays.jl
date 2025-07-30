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
        FarfieldSource is a struct for patterns imported fron CST Studio
"""
struct FarfieldSource <: ElementPattern 
    theta::Vector{Float64}
    phi::Vector{Float64}
    Eθ::Matrix{ComplexF64}
    Eϕ::Matrix{ComplexF64}
end


function FarfieldSource(filename::String)::FarfieldSource
    data = parse_ffs(filename)

    theta = unique(data["theta"])
    phi = unique(data["phi"])
    Nθ = length(theta)
    Nϕ = length(phi)
    Eθ = reshape(data["E_Theta"], length(theta), length(phi))
    Eϕ = reshape(data["E_Phi"], length(theta), length(phi))

    el_patt = FarfieldSource(
        theta,
        phi,
        Eθ,
        Eϕ
    )
    return el_patt
end

function FarfieldSource(data::Dict)::FarfieldSource
    theta = unique(data["theta"])
    phi = unique(data["phi"])
    Nθ = length(theta)
    Nϕ = length(phi)
    Eθ = reshape(data["E_Theta"], length(theta), length(phi))
    Eϕ = reshape(data["E_Phi"], length(theta), length(phi))

    el_patt = FarfieldSource(
        theta,
        phi,
        Eθ,
        Eϕ
    )
    return el_patt
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

function element_gain(p::FarfieldSource, θ::Number, ϕ::Number)
    i = findfirst(p.theta .== θ)
    j = findfirst(p.phi .== ϕ)

    return p.Eθ[i,j], p.Eϕ[i,j]
end

##
##  misc
##
function interpolate_ffs(p::FarfieldSource,theta_grid,phi_grid)::FarfieldSource
    if p.theta == theta_grid && p.phi == phi_grid
        return p
    end
    theta = p.theta
    phi = p.phi
    Eθ = p.Eθ
    Eϕ = p.Eϕ
    Eθ_new = zeros(ComplexF64,length(theta_grid),length(phi_grid))
    Eϕ_new = zeros(ComplexF64,length(theta_grid),length(phi_grid))
    for (i,θ) in enumerate(theta_grid)
        for (j,ϕ) in enumerate(phi_grid)
            theta_idx = argmin(abs.(theta .- θ))
            if theta[theta_idx] - θ < 0
                theta_idx1 = theta_idx
                theta_idx2 = theta_idx+1
            elseif theta[theta_idx] - θ > 0
                theta_idx1 = theta_idx-1
                theta_idx2 = theta_idx
            else
                theta_idx1 = theta_idx2 = theta_idx
            end
    
            phi_idx = argmin(abs.(phi .- ϕ))
            if phi[phi_idx] - ϕ < 0
                phi_idx1 = phi_idx
                phi_idx2 = phi_idx+1
            elseif phi[phi_idx] - ϕ > 0
                phi_idx1 = phi_idx-1
                phi_idx2 = phi_idx
            else
                phi_idx1 = phi_idx2 = phi_idx
            end
            u =  (theta[theta_idx2] - theta[theta_idx1]) == 0 ? 0 : (θ - theta[theta_idx1]) / (theta[theta_idx2] - theta[theta_idx1])
            v = (phi[phi_idx2] - phi[phi_idx1]) == 0 ? 0 : (ϕ - phi[phi_idx1]) / (phi[phi_idx2] - phi[phi_idx1])
            z00 = Eθ[theta_idx1,phi_idx1]
            z01 = Eθ[theta_idx1,phi_idx2]
            z10 = Eθ[theta_idx2,phi_idx1]
            z11 = Eθ[theta_idx2,phi_idx2]
            Eθ_new[i,j] = (z00*(1-u) + z10*u)*(1-v)+
                            (z01*(1-u)+z11*u)*v
            z00 = Eϕ[theta_idx1,phi_idx1]
            z01 = Eϕ[theta_idx1,phi_idx2]
            z10 = Eϕ[theta_idx2,phi_idx1]
            z11 = Eϕ[theta_idx2,phi_idx2]
            Eϕ_new[i,j] = (z00*(1-u) + z10*u)*(1-v)+
                            (z01*(1-u)+z11*u)*v
        end
    end
    return FarfieldSource(
        theta_grid,
        phi_grid,
        Eθ_new,
        Eϕ_new
    )

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
