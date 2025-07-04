struct Pattern
    theta
    phi
    Υθ
    Υϕ
    Υ
end

function Pattern(theta,phi,Υθ,Υϕ)
        Υ = sqrt.(abs.(Υθ) .^2 + abs.(Υϕ).^2)
        return Pattern(theta,phi,Υθ,Υϕ,Υ)
end


function calculate_delays(Array, plane_wave)
    p = Array.coordinates 
    a = plane_wave.a
    return p*a
end

function calculate_array_manifold(Array, plane_wave::PlaneWave)
    p = Array.coordinates
    k = plane_wave.k
    Vk = exp.(-1im .* p*k)
    return Vk
end

function calculate_response(Array, plane_wave::PlaneWave)
    
    w = Array.weights
    Vk = calculate_array_manifold(Array, plane_wave)

    Gθ,Gϕ = element_gain(Array.element_pattern, plane_wave.θ, plane_wave.ϕ)
    F = Gθ * plane_wave.Eθ + Gϕ * plane_wave.Eϕ

    return dot(w, Vk) * F
end

function calculate_pattern(Array,freq)
    theta = 0:0.5:180
    phi = 0:360
    Υθ = zeros(ComplexF64,(length(theta),length(phi)))
    Υϕ = zeros(ComplexF64,(length(theta),length(phi)))
    w = Array.weights
    p = Array.coordinates
    β = 2π * freq / 3e8
    k = [[-β * sind(θ) * cosd(ϕ),-β * sind(θ) * sind(ϕ),-β * cosd(θ)] for θ in theta, ϕ in phi]
    Gθ = [element_gain(Array.element_pattern, θ, ϕ)[1] for θ in theta, ϕ in phi]
    Gϕ = [element_gain(Array.element_pattern, θ, ϕ)[2] for θ in theta, ϕ in phi]
    AF = [sum(w .*cis.(-p*k[i,j])) for i in 1:length(theta),j in 1:length(phi)]
    Υθ = AF .* Gθ
    Υϕ = AF .* Gϕ


    return Pattern(theta,phi,Υθ,Υϕ)
end

##########################
##                      ##
## Pattern statistics   ## 
##                      ##
##########################


function get_power_pattern(pattern::Pattern)
    return abs.(pattern.Υ).^2
end

function get_directivity(pattern::Pattern)
    theta = pattern.theta
    phi = pattern.phi
    P = get_power_pattern(pattern)
    # assuming uniform grid in θ ϕ
    dtheta = diff(theta)[1]/180*π
    dphi = diff(phi)[1]/180*π
    sinθ = sind.(theta)
    W = sinθ .* ones(1,length(phi))
    Prad = sum(P .* W) * dtheta * dphi /4π
    return P ./Prad
end

function get_SLL(pattern::Pattern,threshold_low_dB=-150,threshold_high_db=-10)
    main_beam = maximum(pattern.Υ)
    MB_idx = findall(x -> isapprox(x, MB), pattern.Υ)

    P_dB = 20 * log10.(pattern.Υ ./ MB)

    Mask = falses(size(P_dB))
    visited = falses(size(P_dB)) 

    for start in MB_idx
        if visited[start]
            continue
        end

        queue = [start]

        while !isempty(queue)
            current = popfirst!(queue)

            if visited[current]
                continue
            end

            visited[current] = true
            Mask[current] = true

            current_level = P_dB[current]

            for di in -1:1, dj in -1:1
                neighbor = CartesianIndex(Tuple(current) .+ (di, dj))

                # skip out-of-bounds
                if !checkbounds(Bool, P_dB, neighbor)
                    continue
                end

                neighbor_level = P_dB[neighbor]

                if !visited[neighbor] && neighbor_level > threshold_low_dB
                    if current_level ≤ threshold_high_db 
                        if neighbor_level ≤ current_level   
                            push!(queue, neighbor)
                        end
                    else
                        push!(queue, neighbor)
                    end
                end
            end
        end
    end
    P_dB[Mask] .= -Inf
    return maximum(P_dB)
end

function get_main_beam_direction(pattern::Pattern)

    theta = pattern.theta
    phi = pattern.phi
    MB = maximum(pattern.Υ)
    MB_idx = findall(x -> isapprox(x, MB), pattern.Υ)
    return [(theta[I[1]], phi[I[2]]) for I in MB_idx]

end

function get_radiated_power(pattern::Pattern)
    theta = pattern.theta
    phi = pattern.phi
    P = GetPowerPattern(pattern)
    # assuming uniform grid in θ ϕ
    dtheta = diff(theta)[1]/180*π
    dphi = diff(phi)[1]/180*π
    sinθ = sind.(theta)
    W = sinθ .* ones(1,length(phi))
    Prad = sum(P .* W) * dtheta * dphi /4π
    return Prad
end
