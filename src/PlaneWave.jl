############################
##                        ##
## Plane wave definition  ##
##                        ##
############################

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

