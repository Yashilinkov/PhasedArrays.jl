#######################
##                   ##
##  ULA definition   ##
##                   ##
#######################
"""
    ULA

Uniform Linear Array (ULA) of antenna elements.

# Fields
- `axis::Char`: The axis along which the array is aligned (`'x'`, `'y'`, or `'z'`).
- `N_elements::Int64`: Number of antenna elements in the array.
- `dist::Float64`: Inter-element spacing in meters.
- `coordinates::Matrix{Float64}`: A 3×N matrix of element coordinates in 3D space.
- `weights`: Complex excitation weights for each element (e.g., a vector of `ComplexF64` values).
- `element_pattern::ElementPattern`: Radiation pattern model for the array elements (e.g., isotropic or custom pattern).

# Description
The `ULA` type defines a uniform linear antenna array with equally spaced elements aligned along a specified Cartesian axis. The array supports arbitrary excitation weights and element radiation patterns, and can be used in array factor calculations, beamforming, and electromagnetic simulations.
"""
mutable struct ULA
    axis::Char
    N_elements::Int64
    dist::Float64
    coordinates::Matrix{Float64}
    weights
    element_pattern::ElementPattern
end


function ULA(axis::Char, N::Int, dist::Float64, weights=nothing, pattern=IsotropicPattern('θ'))
    if !(axis in ('x', 'y', 'z'))
        throw(ArgumentError("axis must be 'x', 'y', or 'z'"))
    end
    
    if axis == 'x'
        x = collect(0.0:dist:(N-1)*dist)
        x .-= mean(x)
        y = zeros(Float64,N)
        z = zeros(Float64,N)
    elseif axis == 'y'
        x = zeros(Float64,N)
        y = collect(0.0:dist:(N-1)*dist)
        y .-= mean(y)
        z = zeros(Float64,N)

    elseif axis == 'z'
        x = zeros(Float64,N)
        y = zeros(Float64,N)
        z = collect(0.0:dist:(N-1)*dist)
        z .-= mean(z)

    end
    coords = hcat(x,y,z)

    if weights === nothing
        weights = ones(N)./N
    end

    return ULA(axis, N, dist, coords, weights, pattern)
end

#######################
##                   ##
##  URA definition   ##
##                   ##
#######################


mutable struct URA
    plane::String
    N1::Int64
    N2::Int64
    d1::Float64
    d2::Float64
    coordinates::Matrix{Float64}
    weights
    element_pattern::ElementPattern
end

#
# URA constructor
#

function URA(plane::String, N1, N2, d1, d2,weights=nothing, pattern=IsotropicPattern('θ'))
    if !(plane in ("xy","xz","yz"))
        throw(ArgumentError("plane must be xy, xz or yz"))
    end
    
    idx1 = collect(0:N1-1) .- (N1-1)/2
    idx2 = collect(0:N2-1) .- (N2-1)/2

    x = Float64[]
    y = Float64[]
    z = Float64[]

    for i2 in idx2, i1 in idx1
        if plane == "xy"
            push!(x, i1 * d1)
            push!(y, i2 * d2)
            push!(z, 0.0)
        elseif plane == "xz"
            push!(x, i1 * d1)
            push!(y, 0.0)
            push!(z, i2 * d2)
        elseif plane == "yz"
            push!(x, 0.0)
            push!(y, i1 * d1)
            push!(z, i2 * d2)
        end
    end

    coords = hcat(x, y, z)  # 3×N
    N = N1 * N2
    if weights === nothing
        weights = ones(N) ./ N
    end
    return URA(plane, N1, N2, d1,d2, coords, weights,pattern)
end


#######################
##                   ##
##    Ring Array     ##
##    definition     ##
##                   ##
#######################


mutable struct RingArray
    plane::String
    radius::Float64
    N_elements::Int64
    coordinates::Matrix{Float64}
    weights
    element_pattern::ElementPattern
end

function RingArray(plane::String,
            radius,N_elements,
            weights=nothing,
            pattern=IsotropicPattern('θ'))

    if !(plane in ("xy","xz","yz"))
        throw(ArgumentError("plane must be xy, xz or yz"))
    end
    Δϕ = 360/(N_elements)
    ϕ = 0:Δϕ:360-Δϕ
    if plane == "xy"
        x = radius.*cosd.(ϕ)
        y = radius.*sind.(ϕ)
        z = zeros(N_elements)
    elseif plane == "xz"
        x = radius.*cosd.(ϕ)
        y = zeros(N_elements)
        z = radius.*sind.(ϕ)
    elseif plane == "yz"
        x = zeros(N_elements)
        y = radius.*cosd.(ϕ)
        z = radius.*sind.(ϕ)
    end
    coords = hcat(x, y, z)
    if weights === nothing
        weights = ones(N) ./ N
    end
    return RingArray(plane,radius,N_elements,coords,weights,pattern)
end

#######################
##                   ##
##  Circular Array   ##
##    definition     ##
##                   ##
#######################


@enum LatticeType begin
    Rectangular
    
    
end

mutable struct CircularArray
    plane::String
    R::Float64
    N_elements::Int64
    lattice::LatticeType
    coordinates::Matrix{Float64}
    weights::Vector{ComplexF64}
    element_pattern::ElementPattern
end

function CircularArray(plane::String, R::Float64, lattice::LatticeType, args...; kwargs...)
    if lattice == Rectangular
        return CircularArray(plane, R, Val(:rectangular), args...; kwargs...)
    else
        error("Lattice type $lattice not supported yet")
    end
end



function CircularArray(plane::String, 
            R::Float64, 
            ::Val{:rectangular}, 
            d1::Float64, 
            d2::Float64,
            weights=nothing,
            element_pattern=IsotropicPattern('θ'))
    N1 = 2R÷d1
    N2 = 2R÷d2
    if plane == "xy"
        x = collect( 0:(N1))*d1
        y = collect( 0:(N2))*d2
        z = [0.0]
        x .-= mean(x)
        y .-= mean(y)
    elseif plane == "xz"
        x = collect( 0:(N1))*d1
        y = [0.0]
        z = collect( 0:(N2))*d2
        x .-= mean(x)
        z .-= mean(z)
    elseif plane == "yz"
        x = [0.0]
        y = collect( 0:(N1))*d1
        z = collect( 0:(N2))*d2
        y .-= mean(y)
        z .-= mean(z)
    end
    r = [hypot(x,y,z)/R for x=x,y=y,z=z]
    mask = r .≤ 1 + eps()
    coords = [[x[i],y[j],z[k]] for i in eachindex(x) , j in eachindex(y), k in eachindex(z) if mask[i,j,k]]
    N_elements = sum(mask)

    coords = Matrix(transpose(hcat(coords...)))

    if isnothing(weights)
        weights = ones(ComplexF64,N_elements) ./ N_elements
    end
    return CircularArray(plane, 
            R, 
            N_elements,
            Rectangular, 
            coords, 
            weights, 
            element_pattern)
end
