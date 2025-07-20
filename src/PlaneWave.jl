############################
##                        ##
## Plane wave definition  ##
##                        ##
############################
"""
    PlaneWave(θ, ϕ, frequency; Eθ=1.0 + 0im, Eϕ=0.0 + 0im) -> PlaneWave

Construct a `PlaneWave` object representing an incident electromagnetic wave.

The wave is defined by its direction `(θ, ϕ)` in spherical coordinates (in degrees),
its frequency (in Hz), and its electric field polarization components in the local θ and ϕ directions.

# Arguments
- `θ::Real`: Elevation angle in degrees (0 = +z, 90 = xy-plane).
- `ϕ::Real`: Azimuth angle in degrees (0 = +x, 90 = +y).
- `frequency::Real`: Frequency of the plane wave in Hz.
- `Eθ::ComplexF64` (optional): Complex amplitude of the electric field in θ-direction. Default = `1.0 + 0im`.
- `Eϕ::ComplexF64` (optional): Complex amplitude of the electric field in ϕ-direction. Default = `0.0 + 0im`.

# Returns
- `PlaneWave` object with fields:
  - `θ`, `ϕ`: Direction of arrival (degrees)
  - `frequency`: Frequency (Hz)
  - `a::Vector{Float64}`: Unit direction vector `[ax, ay, az]` (arrival direction)
  - `k::Vector{Float64}`: Wave vector (rad/m)
  - `Eθ`, `Eϕ`: Complex electric field components

# Example
```
pw = PlaneWave(30, 45, 2.4e9)
```
"""
struct PlaneWave
    θ::Number
    ϕ::Number
    frequency::Number
    a::Vector{Float64}
    k::Vector{Float64}
    Eθ::ComplexF64
    Eϕ::ComplexF64
end

function PlaneWave(θ,ϕ,frequency,Eθ::ComplexF64=1.0+0.0im,Eϕ::ComplexF64=0.0+0.0im)
    a_x = -sind(θ)*cosd(ϕ)
    a_y = -sind(θ)*sind(ϕ)
    a_z = -cosd(θ)
    a = [a_x,a_y,a_z]
    k = 2*π*frequency/3e8 .*a
    return PlaneWave(θ,ϕ,frequency,a,k,Eθ,Eϕ)
end

