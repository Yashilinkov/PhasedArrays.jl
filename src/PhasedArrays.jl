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

# ------- exports ----------
export PlaneWave

end
