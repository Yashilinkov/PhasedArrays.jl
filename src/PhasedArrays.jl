module PhasedArrays

# ------- imports --------- 
using Statistics
using LinearAlgebra
using SpecialFunctions
# using GLMakie
# using Base.Threads
using StaticArrays

# -------- includes -------

include("PlaneWave.jl")
include("ElementPatterns.jl")
include("Arrays.jl")
include("Pattern.jl")
include("Beamforming.jl")

# ------- exports ----------
export PlaneWave
export TabulatedPattern,
    IsotropicPattern,
    element_gain,
    compute_dipole_pattern
export ULA, URA
export Pattern, 
    calculate_delays,
    calculate_array_manifold,
    calculate_response,
    calculate_pattern,
    get_power_pattern,
    get_directivity,
    get_SLL
export HansenWoodyardArray
export apply_taper_2D!,
    linear_hamming_window,
    linear_hann_window,
    linear_bartlett_hann_window,
    linear_blackman_window,
    linear_blackman_harris_window,
    linear_kaiser_window,
    linear_tukey_window,
    linear_gaussian_window,
    linear_dolph_cheb_window,
    linear_taylor_window,
    apply_radial_taper!,
    radial_hann_window,
    radial_hamming_window,
    Radial_bartlett_hann_window,
    radial_blackman_window,
    radial_blackman_harris_window,
    radial_kaiser_window,
    radial_tukey_window
export woodward_sampling,
    least_squares_synthesis,
    flattop,
    raised_cosine,
    sinc_pattern,
    gaussian_pattern


end
