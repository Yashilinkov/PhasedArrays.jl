# PhasedArrays.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://Yashilinkov.github.io/PhasedArrays.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Yashilinkov.github.io/PhasedArrays.jl/dev/)
[![Build Status](https://github.com/Yashilinkov/PhasedArrays.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Yashilinkov/PhasedArrays.jl/actions/workflows/CI.yml?query=branch%3Amain)


**PhasedArrays.jl** is a Julia package for simulation and analysis of phased antenna arrays. It supports uniform linear and rectangular arrays, for which it can calculate radiation patterns, also a collection of basic taper functions available. For URA both separable tapers and tapers with radial symmetry are available. for ULA two beamforming schemes provided -- Woodward sampling and least squares synthesis. 

user-defined element patterns (same for all elements) in $\theta,\phi$ coordinates can be defined. 

This is still work in progress.

---
## How to use

assuming Julia is already installed do the following in repl:

```
]
add https://github.com/Yashilinkov/PhasedArrays.jl.git
```

then you can use the package with
```
using PhasedArrays
```
