#######################
##                   ##
##  ULA definition   ##
##                   ##
#######################

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