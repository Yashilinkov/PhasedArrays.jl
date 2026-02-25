"""
    Pattern(theta, phi, Υθ, Υϕ)

Construct a `Pattern` object representing the far-field response of an array.

# Arguments
- `theta::AbstractVector{<:Real}`: Vector of θ angles (degrees).
- `phi::AbstractVector{<:Real}`: Vector of ϕ angles (degrees).
- `Υθ::Matrix{ComplexF64}`: Complex θ-polarized far-field response.
- `Υϕ::Matrix{ComplexF64}`: Complex ϕ-polarized far-field response.

# Returns
A `Pattern` object with total complex field:
    Υ = sqrt.(abs2.(Υθ) + abs2.(Υϕ))
"""

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

"""
    calculate_array_manifold(Array, plane_wave::PlaneWave)

Compute the array manifold vector for a given plane wave.

# Returns
A complex vector of phase shifts for all array elements.
"""
function calculate_array_manifold(Array, plane_wave::PlaneWave)
    p = Array.coordinates
    k = plane_wave.k
    Vk = exp.(-1im .* p*k)
    return Vk
end

"""
    calculate_response(Array, plane_wave::PlaneWave)

Compute complex array response (scalar) to an incoming `PlaneWave`.

Combines:
- Array weights
- Element pattern
- Plane wave polarization and direction

# Returns
Response of the array with given element pattern and weights.
"""
function calculate_response(Array, plane_wave::PlaneWave)
    
    w = Array.weights
    Vk = calculate_array_manifold(Array, plane_wave)

    Gθ,Gϕ = element_gain(Array.element_pattern, plane_wave.θ, plane_wave.ϕ)
    F = Gθ * plane_wave.Eθ + Gϕ * plane_wave.Eϕ

    return dot(w, Vk) * F
end

"""
    calculate_pattern(arr, freq) -> Pattern

Compute the far-field radiation pattern of an antenna array at a given frequency.

This function evaluates the complex far-field response over a spherical grid of angles `(θ, ϕ)`
assuming plane-wave excitation and using the array weights and element pattern. The result includes
the θ- and ϕ-polarized components as well as the combined total field magnitude.

# Arguments
- `arr`: A struct representing the array.
- `freq::Real`: Frequency in Hz at which to compute the pattern.

# Returns
- `Pattern` object containing:
  - `theta`: θ angles in degrees (`0:0.5:180`).
  - `phi`: ϕ angles in degrees (`0:1:360`).
  - `Υθ`: θ-polarized complex far-field.
  - `Υϕ`: ϕ-polarized complex far-field.
  - `Υ`: Total field magnitude √(|Υθ|² + |Υϕ|²)
"""
function calculate_pattern(arr, freq)
    theta = 0:0.5:180
    phi = 0:360
    Nθ, Nϕ = length(theta), length(phi)

    Υθ = Matrix{ComplexF64}(undef, Nθ, Nϕ)
    Υϕ = Matrix{ComplexF64}(undef, Nθ, Nϕ)

    w = arr.weights
    p = arr.coordinates

    N = arr.N_elements
    β = 2π * freq / 3e8

    sinθarr = sind.(theta)
    cosθarr = cosd.(theta)
    sinϕarr = sind.(phi)
    cosϕarr = cosd.(phi)

    @inbounds for i in 1:Nθ
        sinθ = sinθarr[i]
        cosθ = cosθarr[i]
        @inbounds for j in 1:Nϕ
            sinϕ = sinϕarr[j]
            cosϕ = cosϕarr[j]

            kx = -β * sinθ * cosϕ
            ky = -β * sinθ * sinϕ
            kz = -β * cosθ

            acc = 0.0 + 0.0im
            @inbounds for n in 1:N
                dotprod = p[n,1]*kx + p[n,2]*ky + p[n,3]*kz
                acc += w[n] * cis(-dotprod)
            end

            Υθ[i,j] = acc * element_gain(arr.element_pattern, theta[i], phi[j])[1] 
            Υϕ[i,j] = acc * element_gain(arr.element_pattern, theta[i], phi[j])[2] 
        end
    end

    return Pattern(theta, phi, Υθ, Υϕ)
end

function calculate_pattern(arr::CustomPhasedArray)
    theta = 0:0.5:180
    phi = 0:360
    interpolate_custom_array!(arr,theta,phi)
    Nθ, Nϕ = length(theta), length(phi)

    Υθ = zeros(ComplexF64,Nθ,Nϕ)
    Υϕ = zeros(ComplexF64,Nθ,Nϕ)

    w = arr.weights
    # phase_centres = arr.phase_centres

    N = arr.N_elements

    @inbounds for n in 1:N
        Υθ .+= arr.element_patterns[n].Eθ .* w[n]
        Υϕ .+= arr.element_patterns[n].Eϕ .* w[n]
                
    end

    return Pattern(theta, phi, Υθ, Υϕ)
end

function interpolate_custom_array!(arr::CustomPhasedArray, theta_grid, phi_grid)
    # materialize to vectors (Ranges are fine, but we’ll use values a lot)
    θg = collect(theta_grid)
    ϕg = collect(phi_grid)

    # If user gives 0:360, drop the last point (duplicate of 0)
    if !isempty(ϕg) && isapprox(ϕg[end] - ϕg[1], 360; atol=1e-12, rtol=0)
        ϕg = ϕg[1:end-1]
    end

    # Normalize phi grid to [0,360) and keep it sorted (important for argmin neighbor logic)
    ϕg = mod.(ϕg,360)
    sort!(ϕg)

    # Interpolate each element’s pattern to θg, ϕg
    @inbounds for n in 1:arr.N_elements
        arr.element_patterns[n] = interpolate_ffs(arr.element_patterns[n], θg, ϕg)
    end

    return arr
end



##########################
##                      ##
## Pattern statistics   ## 
##                      ##
##########################


function get_power_pattern(pattern::Pattern)
    return abs2.(pattern.Υ)
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
    MB_idx = findall(x -> isapprox(x, main_beam), pattern.Υ)

    P_dB = 20 * log10.(pattern.Υ ./ main_beam)

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
    P = get_power_pattern(pattern)
    # assuming uniform grid in θ ϕ
    dtheta = diff(theta)[1]/180*π
    dphi = diff(phi)[1]/180*π
    sinθ = sind.(theta)
    W = sinθ .* ones(1,length(phi))
    Prad = sum(P .* W) * dtheta * dphi /4π
    return Prad
end
