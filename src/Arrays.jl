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
The `ULA` type defines a uniform linear antenna array with equally spaced elements aligned along a specified Cartesian axis. 
The array supports arbitrary excitation weights and element radiation patterns.

Example use:
```@example
using PhasedArrays
f = 3e9
λ = 3e8/f
d = 0.5*λ
N = 5
ula = ULA('x',N,d)
```
"""
mutable struct ULA
    axis::Char
    N_elements::Int64
    dist::Float64
    coordinates::Matrix{Float64}
    weights::Vector{ComplexF64}
    element_pattern::ElementPattern
end


"""
    ULA(axis::Char, N::Int, dist::Float64; weights=nothing, pattern=IsotropicPattern('θ'))

Construct a Uniform Linear Array along the given `axis` (`'x'`, `'y'`, or `'z'`) with `N` elements spaced `dist` meters apart.

Weights and element pattern are optional; defaults are uniform weighting and an isotropic θ-directed pattern.
"""
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

"""
    URA

Uniform Rectangular Array (URA) of antenna elements.

# Fields
- `plane::String`: The Cartesian plane in which the array lies. Must be `"xy"`, `"xz"`, or `"yz"`.
- `N1::Int64`: Number of elements along the first axis of the plane.
- `N2::Int64`: Number of elements along the second axis of the plane.
- `d1::Float64`: Element spacing along the first axis (in meters).
- `d2::Float64`: Element spacing along the second axis (in meters).
- `weights`: Complex excitation weights for each element (e.g., a vector of `ComplexF64` values).
- `element_pattern::ElementPattern`: Radiation pattern of the element.
# Description
The `URA` type represents a uniform rectangular antenna array, constructed in one of the three primary Cartesian planes.
It supports arbitrary element weights and customizable radiation patterns.

Example use:
```@example
using PhasedArrays
f = 3e9
λ = 3e8/f
dx = 0.5*λ
dy = 0.25*λ
N = 5
ura = URA("xy",N,N,dx,dy)
```
"""
mutable struct URA
    plane::String
    N1::Int64
    N2::Int64
    d1::Float64
    d2::Float64
    coordinates::Matrix{Float64}
    weights::Vector{ComplexF64}
    element_pattern::ElementPattern
end

#
# URA constructor
#

"""
    URA(plane::String, N1, N2, d1, d2; weights=nothing, pattern=IsotropicPattern('θ'))

Construct a Uniform Rectangular Array (URA) with `N1 × N2` elements spaced by `d1` and `d2` meters along the specified `plane`.

# Arguments
- `plane`: The plane in which the array is placed. Must be `"xy"`, `"xz"`, or `"yz"`.
- `N1`: Number of elements along the first axis of the chosen plane.
- `N2`: Number of elements along the second axis of the chosen plane.
- `d1`: Spacing along the first axis (in meters).
- `d2`: Spacing along the second axis (in meters).
- `weights` (optional): Excitation weights; defaults to uniform weights (`1/N`).
- `pattern` (optional): Element radiation pattern; defaults to isotropic pattern in the θ-direction.

# Returns
A `URA` instance with 3D element coordinates centered around the origin.

# Example
```julia
using PhasedArrays
f = 10e9
λ = 3e8/f
ura = URA("xy", 4, 4, 0.5λ, 0.5λ)
```
"""
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

"""
    RingArray

Circular (ring) array of antenna elements in a specified Cartesian plane.

# Fields
- `plane::String`: The plane in which the ring lies. Must be `"xy"`, `"xz"`, or `"yz"`.
- `radius::Float64`: Radius of the ring in meters.
- `N_elements::Int64`: Number of equally spaced elements along the ring.
- `coordinates::Matrix{Float64}`: A 3×N matrix of element coordinates in 3D space.
- `weights`: Complex excitation weights for each element (e.g., a vector of `ComplexF64` values).
- `element_pattern::ElementPattern`: Radiation pattern model for the array elements (e.g., isotropic or directional).

# Description
The `RingArray` type defines a circular array with `N_elements` evenly spaced antennas placed on a ring of radius `radius` in the specified plane.
It supports arbitrary excitation weights and customizable element radiation patterns.
"""
mutable struct RingArray
    plane::String
    radius::Float64
    N_elements::Int64
    coordinates::Matrix{Float64}
    weights::Vector{ComplexF64}
    element_pattern::ElementPattern
end

"""
    RingArray(plane::String, radius, N_elements; weights=nothing, pattern=IsotropicPattern('θ'))

Construct a circular antenna array with `N_elements` elements equally spaced on a ring of radius `radius` in the specified `plane`.

# Arguments
- `plane`: The Cartesian plane in which the ring lies (`"xy"`, `"xz"`, or `"yz"`).
- `radius`: Radius of the ring in meters.
- `N_elements`: Number of antenna elements placed uniformly on the ring.
- `weights` (optional): Excitation weights for each element; defaults to uniform weights (`1/N`).
- `pattern` (optional): Element radiation pattern; defaults to an isotropic θ-directed pattern.

# Returns
A `RingArray` instance with 3D element coordinates arranged around the origin.

# Example
```julia
using PhasedArrays
f = 5e9
λ = 3e8/f
ra = RingArray("xy", λ, 8)
```
"""
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

"""
    CircularArray

Circular aperture array with support for various lattice types.

# Fields
- `plane::String`: The Cartesian plane in which the circular aperture lies (`"xy"`, `"xz"`, or `"yz"`).
- `R::Float64`: Radius of the circular aperture in meters.
- `N_elements::Int64`: Number of antenna elements in the array.
- `lattice::LatticeType`: Lattice geometry used for placing elements (e.g., `Rectangular`).
- `coordinates::Matrix{Float64}`: A 3×N matrix of element coordinates in 3D space.
- `weights::Vector{ComplexF64}`: Excitation weights applied to each element.
- `element_pattern::ElementPattern`: Radiation pattern model of each element.

# Description
The `CircularArray` type represents a 2D array of antenna elements placed within a circular aperture of radius `R`, lying in the specified Cartesian plane.
The layout of elements is controlled by the chosen `lattice` geometry (currently, only `Rectangular` is supported).
This array type supports arbitrary excitation weights and element radiation patterns.
"""
mutable struct CircularArray
    plane::String
    R::Float64
    N_elements::Int64
    lattice::LatticeType
    coordinates::Matrix{Float64}
    weights::Vector{ComplexF64}
    element_pattern::ElementPattern
end

"""
    CircularArray(plane::String, R::Float64, lattice::LatticeType, args...; kwargs...)

Construct a circular array of elements within a radius `R` in the specified Cartesian `plane` using the specified `lattice` type.

# Arguments
- `plane`: Plane in which the array lies (`"xy"`, `"xz"`, or `"yz"`).
- `R`: Radius of the circular aperture (in meters).
- `lattice`: Lattice type used for element placement (currently only `Rectangular` is supported).
- `args...`, `kwargs...`: Passed to the specific lattice constructor (e.g., inter-element spacing, weights, element pattern).

# Returns
A `CircularArray` with coordinates generated based on the chosen lattice type.

# Example
```julia
f = 3e9
λ = 3e8 / f
array = CircularArray("xy", 5λ, Rectangular, 0.5λ, 0.5λ)
```
"""
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
