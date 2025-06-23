##########################
##                      ##
## Tapers               ##
##                      ##
##########################

function HansenWoodyardArray(axis::Char, N::Int, freq)
    λ = 3e8/freq
    d = λ/2 * (1-1/N)
    array = ULA('x',N,d)
    pw = PlaneWave(90.0,0.0,freq)

    w = calculate_array_manifold(array,pw)
    array.axis = axis
    array.weights = w
    return array
end


##                      ##
## Separable tapers     ##
##  (1d tapers)         ##


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
    # Construct 2D taper by outer product (N2 x N1)
    weights2D = w1 * transpose(w2)
    # Reshape weights vector into matrix (N2 rows, N1 cols)
    weightsMatrix = reshape(Array.weights, (Array.N1, Array.N2))

    weightsMatrix .*= weights2D

    # Update weights vector
    Array.weights .= vec(weightsMatrix) ./ sum(abs, weightsMatrix)
    return nothing
end



function linear_hamming_window(N)
    w = 0.54 .- 0.46 .* cos.(2π .* (0:N-1) ./ (N-1))
    return w./ maximum(abs.(w))
end

function linear_hann_window(N)
    w = 0.5 .- 0.5 .* cos.(2π .* (0:N-1) ./ (N-1))
    return w./ maximum(abs.(w))
end

function linear_bartlett_hann_window(N)

    w = zeros(Float64, N)
    for n in 0:N-1
        x = n / (N-1) - 0.5
        w[n+1] = 0.62 - 0.48 * abs(x) + 0.38 * cos(2π * x)
    end
    return w./ maximum(abs.(w))
end

function linear_blackman_window(N)
    n = 0:N-1
    w = 0.42 .- 0.5*cos.(2π*n/(N-1)) .+ 0.08*cos.(4π*n/(N-1))
    return w./ maximum(abs.(w))
end

function linear_blackman_harris_window(N)
    n = 0:N-1
    a0, a1, a2, a3 = 0.35875, 0.48829, 0.14128, 0.01168
    w = a0 .- a1*cos.(2π*n/(N-1)) .+ a2*cos.(4π*n/(N-1)) .- a3*cos.(6π*n/(N-1))
    return w./ maximum(abs.(w))
end

function linear_kaiser_window(N, β)
    n = 0:N-1
    α = (N-1)/2
    w = besseli.(0, β .* sqrt.(1 .- ((n .- α) ./ α).^2)) ./ besseli(0, β)
    return w./ maximum(abs.(w))
end

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

function linear_gaussian_window(N, σ)
    n = 0:N-1
    a = (N - 1) / 2
    w = exp.(-0.5 * ((n .- a) ./ (σ * a)) .^ 2)
    return w./ maximum(abs.(w))
end

function linear_dolph_cheb_window(N,R,d, λ)
    # N - number of elements
    # R - SLL in dB (for example -20)
    # d interelement distance (m)
    # λ wavelength (m)
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

function linear_taylor_window(N,R,n̄,d,λ)
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

function radial_hann_window(r)
    w = 0.5 * (1 - cos(π * (1 - r)))
    return w
end

function radial_hamming_window(r)
    return 0.54 - 0.46 * cos(2π * (0.5r+0.5))
end

function Radial_bartlett_hann_window(r)
    return 0.62 - 0.48 * abs(0.5*r) + 0.38 * cos(π * r)
end

function radial_blackman_window(r)
    return 0.42 - 0.5*cos(2π*(0.5r+0.5)) + 0.08*cos(4π*(0.5r+0.5))
end

function radial_blackman_harris_window(r)
    a₀, a₁, a₂, a₃ = 0.35875, 0.48829, 0.14128, 0.01168
    r̄ = 0.5*r-0.5
   return a₀ - a₁*cos(2π*r̄)+a₂*cos(4π*r̄) - a₃*cos(6π*r̄)
end

function radial_kaiser_window(r,β)
    if r ≤ 1
        return besseli(0,β*sqrt(1-r^2))/besseli(0, β)
    else
        return 0
    end
end


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

function flattop(HalfWidth)
    # half width in degrees
    return u -> abs(u) ≤ abs(cosd(90-HalfWidth)) ? 1.0 : 0.0
end

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

# usage:
function sinc_pattern(scale)
    return u -> sinc(scale * u)
end


function gaussian_pattern(σ)
    return u -> exp(-(u^2) / (2 * σ^2))
end

