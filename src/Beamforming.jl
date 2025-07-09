##########################
##                      ##
## Tapers               ##
##                      ##
##########################

"""
    HansenWoodyardArray(axis::Char, N::Int, freq::Float64) -> ULA

Constructs a Hansen-Woodyard endfire linear array along the specified axis.

# Arguments
- `axis`: Array orientation; must be one of `'x'`, `'y'`, or `'z'`.
- `N`: Number of antenna elements.
- `freq`: Operating frequency in Hz.

# Returns
A `ULA` object with Hansen-Woodyard spacing and weights that steer the beam toward the endfire direction.

# Notes
- Spacing is set to `λ/2 * (1 - 1/N)` for enhanced endfire directivity.
- The beam is steered to:
  - θ=90°, ϕ=0° for `'x'`-axis
  - θ=90°, ϕ=90° for `'y'`-axis
  - θ=0° for `'z'`-axis

# Example
```julia
array = HansenWoodyardArray('y', 10, 3e9)
```
"""
function HansenWoodyardArray(axis::Char, N::Int, freq)
    if !(axis in ['x', 'y', 'z'])
        throw(ArgumentError("axis must be one of 'x', 'y', or 'z'"))
    end
    
    λ = 3e8/freq
    d = λ/2 * (1-1/N)
    if axis == 'x'
        array = ULA('x',N,d)
        pw = PlaneWave(90.0,0.0,freq)
    elseif axis == 'y'
        array = ULA('y',N,d)
        pw = PlaneWave(90.0,90.0,freq)
    elseif axis == 'z'
        array = ULA('z',N,d)
        pw = PlaneWave(0.0,0.0,freq)
    end
    w = calculate_array_manifold(array,pw)
    array.weights = w
    return array
end


##                      ##
## Separable tapers     ##
##  (1d tapers)         ##

"""
    apply_taper_2D!(array::URA, taper::Function, axis::Symbol)

Applies a 2D tapering function (with separable kernels) to the weights of a uniform rectangular array (`URA`).

# Arguments
- `array`: A `URA` instance whose weights will be modified in-place.
- `taper`: A function that takes an integer `N` and returns a vector of length `N` representing the taper weights along one dimension.
- `axis`: Specifies which dimensions to apply the taper on.  
    - `:dim1` applies taper along the first axis.  
    - `:dim2` applies taper along the second axis.  
    - `:both` applies same taper on both axes.

# Behavior
The function constructs a 2D taper weight matrix by taking the outer product of the taper vectors along each dimension.  
It multiplies the existing weights element-wise by this taper matrix and then normalizes the resulting weights so that the sum of their absolute values is 1.

# Notes
- The `weights` field of `array` is updated in-place.
- The function assumes the weights vector can be reshaped into `(N1, N2)`, where `N1` and `N2` are the array dimensions.

# Example
```julia
tukey = N -> linear_tukey_window(N, 0.5)  # Example taper function
apply_taper_2D!(ura, tukey, :both)
```
"""
function apply_taper_2D!(Array::URA,Taper,axis)
    # Taper is a function that accepts N (length) and returns a vector window of length N
    # axis can be :both, :dim1, or :dim2
    w1 = nothing
    w2 = nothing
    if axis == :dim1 || axis == :both
        w1 = Taper(Array.N1)
    else
        w1 = ones(Array.N1)
    end
    if axis == :dim2 || axis == :both
        w2 = Taper(Array.N2)
    else
        w2 = ones(Array.N2)
    end

    weights2D = w1 * transpose(w2)

    weightsMatrix = reshape(Array.weights, (Array.N1, Array.N2))

    weightsMatrix .*= weights2D


    Array.weights .= vec(weightsMatrix) ./ sum(abs, weightsMatrix)
    return nothing
end


"""
    linear_hamming_window(N)

Returns a normalized Hamming window of length `N`.

Used to taper array weights to reduce sidelobe levels. The window is normalized to have maximum value 1.
"""
function linear_hamming_window(N)
    w = 0.54 .- 0.46 .* cos.(2π .* (0:N-1) ./ (N-1))
    return w./ maximum(abs.(w))
end

"""
    linear_hann_window(N)

Returns a normalized Hann window of length `N`.

A cosine-based taper that reduces sidelobes. Maximum value is normalized to 1.
"""
function linear_hann_window(N)
    w = 0.5 .- 0.5 .* cos.(2π .* (0:N-1) ./ (N-1))
    return w./ maximum(abs.(w))
end

"""
    linear_bartlett_hann_window(N)

Returns a normalized Bartlett-Hann window of length `N`.
Normalized to peak 1.
"""
function linear_bartlett_hann_window(N)

    w = zeros(Float64, N)
    for n in 0:N-1
        x = n / (N-1) - 0.5
        w[n+1] = 0.62 - 0.48 * abs(x) + 0.38 * cos(2π * x)
    end
    return w./ maximum(abs.(w))
end

"""
    linear_blackman_window(N)

Returns a normalized Blackman window of length `N`.

Provides good sidelobe suppression at the cost of wider mainlobe. Normalized to peak 1.
"""
function linear_blackman_window(N)
    n = 0:N-1
    w = 0.42 .- 0.5*cos.(2π*n/(N-1)) .+ 0.08*cos.(4π*n/(N-1))
    return w./ maximum(abs.(w))
end

"""
    linear_blackman_harris_window(N)

Returns a normalized Blackman-Harris window of length `N`.

High sidelobe suppression, typically better than standard Blackman. Normalized to 1.
"""
function linear_blackman_harris_window(N)
    n = 0:N-1
    a0, a1, a2, a3 = 0.35875, 0.48829, 0.14128, 0.01168
    w = a0 .- a1*cos.(2π*n/(N-1)) .+ a2*cos.(4π*n/(N-1)) .- a3*cos.(6π*n/(N-1))
    return w./ maximum(abs.(w))
end

"""
    linear_kaiser_window(N, β)

Returns a normalized Kaiser window of length `N` with shape parameter `β`.

Useful for controlling sidelobe levels via `β`. Normalized to peak 1.
"""
function linear_kaiser_window(N, β)
    n = 0:N-1
    α = (N-1)/2
    w = besseli.(0, β .* sqrt.(1 .- ((n .- α) ./ α).^2)) ./ besseli(0, β)
    return w./ maximum(abs.(w))
end

"""
    linear_tukey_window(N, α)

Returns a normalized Tukey (tapered cosine) window of length `N` with tapering parameter `α` ∈ [0, 1].

- `α = 0` → rectangular window
- `α = 1` → Hann window

Normalized to maximum value 1.
"""
function linear_tukey_window(N, α)
    if α < 0.0 || α > 1.0
        throw(ArgumentError("Tukey window parameter α must be between 0 and 1, got α = $α"))
    end
    n = 0:N-1
    w = zeros(Float64, N)
    for i in eachindex(n)
        t = n[i] / (N-1)
        if t < α/2
            w[i] = 0.5 * (1 + cos(π*(2t/α - 1)))
        elseif t <= 1 - α/2
            w[i] = 1.0
        else
            w[i] = 0.5 * (1 + cos(π*(2t/α - 2/α + 1)))
        end
    end
    return w./ maximum(abs.(w))
end

"""
    linear_gaussian_window(N, σ)

Returns a normalized Gaussian window of length `N` with spread parameter `σ`.

smaller `σ` results in wider mainlobe and lower sidelobes. Normalized to peak 1.
"""
function linear_gaussian_window(N, σ)
    n = 0:N-1
    a = (N - 1) / 2
    w = exp.(-0.5 * ((n .- a) ./ (σ * a)) .^ 2)
    return w./ maximum(abs.(w))
end

"""
    linear_dolph_cheb_window(N, R, d, λ)

Computes Dolph-Chebyshev window of length `N` for a ULA with element spacing `d` and wavelength `λ`.

# Arguments
- `R`: Desired sidelobe level in dB (e.g. -20).
- `d`: Element spacing in meters.
- `λ`: Wavelength in meters.

Returns a real-valued weight vector normalized to unit gain in the endfire direction.
"""
function linear_dolph_cheb_window(N,R)
    λ = 0.01
    d = 0.5λ
    f₀ = 3e8/λ
    ula = ULA('z',N,d)
    R′ = 10^(-R/20)
    x₀ = cosh(1/(N-1)*acosh(R′))
    x_min = x₀*cos(π*d/λ)
    abs(x_min) < 1e-12 && (x_min = 0.0)
    xp = zeros(N-1)
    for p in 1:N-1
        xp[p] = cos( (2p-1)/(N-1) * (π/2) )
    end
    x = xp[x_min .≤ xp .≤ x₀]
    ψ = 2*acos.(x/x₀)
    ψp = vcat(0,ψ, -reverse(ψ) )
    theta = acosd.(ψp*λ/(2π*d))
    Nangles = length(theta)
    V = zeros(ComplexF64,N,length(theta))
    for i in 1:Nangles
        pw = PlaneWave(theta[i],0,f₀)
        v = calculate_array_manifold(ula,pw)
        V[:,i] = v
    end
    e = zeros(length(theta),1)
    e[1] = 1

    w = V' \ e
    return real.(w)
end

"""
    linear_taylor_window(N, R, n̄, d, λ)

Generates a Taylor window with `N` elements, sidelobe level `R` (dB), and `n̄` equal-level sidelobes.

# Arguments
- `R`: Desired sidelobe level in dB (for example -20).
- `n̄`: Number of equal sidelobes + 1. Must be between 1 and N–1.
- `d`: Element spacing in meters.
- `λ`: Wavelength in meters.

Returns a real-valued taper optimized for sidelobe control in antenna arrays.
"""
function linear_taylor_window(N,R,n̄)
    # N - number of elements
    # R - SLL in dB (for example -20)
    # n̄ - number of first equal-ish sidelobes +1
    # d interelement distance (m)
    # λ wavelength (m)
    if n̄ ≥ N
        throw(ArgumentError("n̄ is between 1 and N-1"))
    elseif n̄ ≤ 0
        throw(ArgumentError("n̄ is between 1 and N-1"))
    end
    λ = 0.01
    d = 0.5λ
    ula = ULA('z',N,d)
    R′ = 10^(-R/20)
    x₀ = cosh(1/(N-1)*log(R′+sqrt(R′^2-1)))
    p = collect(1:n̄)
    ψcheb = 2 .* acos.(1/x₀ .* cos.( (2p .- 1).* (π/(2*(N-1))) ) )
    σ = 2π*n̄/(N*ψcheb[n̄])
    ψ′ₙ = σ*ψcheb[1:n̄-1]
    n = collect(n̄:(N-1)/2)
    ψun = 2π/N .* n
    ψtmp = vcat(ψ′ₙ,ψun) .*λ/(2π*d)
    ψ = vcat(0, ψtmp[abs.(ψtmp).<=1], -reverse(ψtmp[abs.(ψtmp).<=1]))
    theta = acosd.(ψ)

    Nangles = length(theta)

    V = zeros(ComplexF64,N,Nangles)
    for i in 1:Nangles
        pw = PlaneWave(theta[i],0,f)
        v = calculate_array_manifold(ula,pw)
        V[:,i] = v
    end
    e = zeros(length(theta),1)
    e[1] = 1

    w = V' \ e

    return real.(w)
end



##                  ##
##  radial tapers   ##
##                  ##

"""
    apply_radial_taper!(array::URA, taper_fun::Function, n::Int)

Applies a radial taper to the weights of a `URA` array in-place.

# Arguments
- `array`: A `URA` (Uniform Rectangular Array) whose weights will be modified.
- `taper_fun`: A function mapping normalized radius ∈ [0, 1] to a taper value.
- `n`: Number of outer radial rings (based on unique distance levels) to suppress fully (i.e., taper to zero).

# Behavior
- Computes distance of each element from array center.
- Applies taper based on normalized radius `r/r_max`, where `r_max` is the distance to the `n`-th inner ring.
- Weights are multiplied by the taper and normalized.

# Example
```julia
apply_radial_taper!(ura, radial_hann_window, 1)
```
"""
function apply_radial_taper!(Array::URA, taper_fun::Function, n)
    # n is the number of rings "removed" from the array, 
    # starting with corner elements
    # n ≥ 0
    coords = eachrow(Array.coordinates)
    N = length(coords)

    # Compute radial distance from array center
    r = zeros(Float64, N)
    for (p_ind,p) in enumerate((coords))
        r[p_ind] = hypot(p...)
    end

    r_sorted = sort(unique(r),rev=true)
    r_max = r_sorted[n+1]
    w = taper_fun.(r ./ r_max)

    mask = r .> r_max
    w[mask] .*= 0

    Array.weights .*= w
    Array.weights ./= sum(abs.(Array.weights))
    return nothing
end

function apply_radial_taper!(Array::CircularArray,taper_fun::Function)
    R = Array.R
    coords = eachrow(Array.coordinates)
    N = Array.N_elements
    w = zeros(ComplexF64,N)
    for (p_ind,p) in enumerate(coords)
        r = hypot(p...)/R
        w[p_ind] = taper_fun(r)
    end
    w ./= sum(abs.(w))
    Array.weights = w
    return nothing
end

"""
    radial_hann_window(r::Float64)

Returns the value of a Hann window evaluated at normalized radius `r ∈ [0, 1]`.

Used in radial tapers for circular or rectangular arrays.
"""
function radial_hann_window(r)
    w = 0.5 * (1 - cos(π * (1 - r)))
    return w
end
"""
    radial_hamming_window(r::Float64)

Returns the value of a Hamming window at normalized radius `r`.

Commonly used for moderate sidelobe suppression in array tapering.
"""
function radial_hamming_window(r)
    return 0.54 - 0.46 * cos(2π * (0.5r+0.5))
end

"""
    Radial_bartlett_hann_window(r::Float64)

Returns Bartlett-Hann window evaluated at normalized radius `r`.

Hybrid taper with smoother roll-off and moderate sidelobe control.
"""
function Radial_bartlett_hann_window(r)
    return 0.62 - 0.48 * abs(0.5*r) + 0.38 * cos(π * r)
end

"""
    radial_blackman_window(r::Float64)

Blackman window evaluated at normalized radius `r`.

Provides stronger sidelobe suppression than Hann or Hamming.
"""
function radial_blackman_window(r)
    return 0.42 - 0.5*cos(2π*(0.5r+0.5)) + 0.08*cos(4π*(0.5r+0.5))
end

"""
    radial_blackman_harris_window(r::Float64)

Returns Blackman-Harris window at normalized radius `r`.

"""
function radial_blackman_harris_window(r)
    a₀, a₁, a₂, a₃ = 0.35875, 0.48829, 0.14128, 0.01168
    r̄ = 0.5*r-0.5
   return a₀ - a₁*cos(2π*r̄)+a₂*cos(4π*r̄) - a₃*cos(6π*r̄)
end

"""
    radial_kaiser_window(r::Float64, β::Float64)

Returns Kaiser window value at normalized radius `r`, with shape parameter `β`.

"""
function radial_kaiser_window(r,β)
    if r ≤ 1
        return besseli(0,β*sqrt(1-r^2))/besseli(0, β)
    else
        return 0
    end
end

"""
    radial_tukey_window(r::Float64, α::Float64)

Returns Tukey (tapered cosine) window value at normalized radius `r`, with tapering parameter `α ∈ [0, 1]`.

- `α = 0` gives rectangular window
- `α = 1` gives Hann window

"""
function radial_tukey_window(r,α)
    if α < 0.0 || α > 1.0
    throw(ArgumentError("Tukey window parameter α must be between 0 and 1, got α = $α"))
    end
    r̄ = 0.5r+0.5
    if r̄ ≤ α/2
        return 0.5 * (1+cos(π*(2(r̄)/α -1 )))
    elseif r̄ ≤ 1-α/2
        return 1
    else
        return 0.5*(1+cos(π*(2*(r̄)/α-2/α+1)))
    end
end


##########################
##                      ##  
##  Beamforming         ##
##                      ##
##########################

"""
    woodward_sampling(N, d, λ, goalPattern) -> Vector{ComplexF64}

Computes array weights using Woodward (spatial) sampling method.

# Arguments
- `N`: Number of elements.
- `d`: Element spacing (in meters).
- `λ`: Wavelength (in meters).
- `goalPattern`: A function `u -> desired pattern`, where `u = cos(θ)` is the directional cosine.

# Returns
Complex weight vector that approximates the desired pattern over sampled points.

"""
function woodward_sampling(N,d,λ,goalPattern)
    # sampling points in directional cosines
    m = collect(0:N-1)
    uₘ= λ/d * 1/N.*(m .- (N-1)/2 )
    k = 2π/λ
    # desired pattern
    Bdu = [goalPattern(u) for u in uₘ]
    # Compute weights via inverse spatial Fourier transform
    w = ComplexF64[]
    w = zeros(ComplexF64,N)
    for n in 1:N
        w[n] = 1/N * sum( Bdu .* exp.(-1im *k * d * uₘ .* ((n - 1) - (N-1)/2)  ))
    end

    return w
    
end

"""
    least_squares_synthesis(N, d, λ, Bd; θrange=0.0:0.1:180.0) -> Vector{ComplexF64}

Computes array weights using least-squares pattern synthesis.

# Arguments
- `N`: Number of elements.
- `d`: Element spacing.
- `λ`: Wavelength.
- `Bd`: Desired pattern as a function of `u = cos(θ)`.
- `θrange`: Angle range in degrees over which the pattern is matched (default: `0.0:0.1:180.0`).

# Returns
Weight vector that minimizes the squared error between the desired and actual beam pattern.

"""
function least_squares_synthesis(N,d,λ,Bd::Function; θrange=0.0:0.1:180.0)
    f = 3e8/λ
    k = 2π/λ

    A = zeros(ComplexF64, N, N)
    I = zeros(ComplexF64, N)
    ula_tmp = ULA('z',N,d)
    for θ in θrange
        pw = PlaneWave(θ, 0.0, f)
        v = CalculateArrayManifold(ula_tmp, pw)
        A .+= v * v' * sind(θ)
        I .+= v * Bd(cosd(θ)) * sind(θ)
    end

    w = A \ I 
    return w
end

#
# example goalPatterns
#
"""
    flattop(HalfWidth)

Returns a function `u -> 1 or 0`, representing a flat-top beam pattern of given half-width (in degrees).

# Arguments
- `HalfWidth`: Angular half-width of the flat region in degrees.

# Returns
Function defined on directional cosine `u = cos(θ)`, with `1` inside the main lobe and `0` outside.
"""
function flattop(HalfWidth)
    # half width in degrees
    return u -> abs(u) ≤ abs(cosd(90-HalfWidth)) ? 1.0 : 0.0
end

"""
    raised_cosine(HalfWidth, rolloff)

Returns a function representing a raised-cosine-shaped beam pattern in the directional cosine domain.

# Arguments
- `HalfWidth`: Half-width of the main lobe (in degrees).
- `rolloff`: Transition width (in directional cosine units).

# Returns
A smooth beam pattern tapering from 1 to 0 around the edge of the flat region.
"""
function raised_cosine(HalfWidth, rolloff)
    u₀ = abs(cosd(90-HalfWidth))
    return u -> begin
        absu = abs(u)
        if absu ≤ u₀ - rolloff
            return 1.0
        elseif absu ≤ u₀
            return 0.5 * (1 + cos(π * (absu - (u₀ - rolloff)) / rolloff))
        else
            return 0.0
        end
    end
end

"""
    sinc_pattern(scale)

Returns a sinc-based beam shaping function.

# Arguments
- `scale`: Controls the lobe spacing and width in the directional cosine domain.

# Returns
A function `u -> sinc(scale * u)`, centered at broadside.
"""
function sinc_pattern(scale)
    return u -> sinc(scale * u)
end

"""
    gaussian_pattern(σ)

Returns a Gaussian-shaped pattern in the directional cosine domain.

# Arguments
- `σ`: Standard deviation controlling the width of the main lobe.

# Returns
A function `u -> exp(-u^2 / (2σ^2))` for smooth, bell-shaped beam shaping.
"""
function gaussian_pattern(σ)
    return u -> exp(-(u^2) / (2 * σ^2))
end



##########################
##                      ##  
##  Misc Beamforming    ##
##      related         ##
##                      ##
##########################

function truncate_phase!(arr,bits::Int)
    if bits == 0
        return nothing
    else
        Δϕ = 2π/2^bits
        w = arr.weights
        phase = round.(angle.(w)/Δϕ).*Δϕ
        arr.weights .= abs.(w).*cis.(phase)
        return nothing   
    end 
end  