# RandomKrylov

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tmigot.github.io/RandomKrylov.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tmigot.github.io/RandomKrylov.jl/dev)
[![Build Status](https://github.com/tmigot/RandomKrylov.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tmigot/RandomKrylov.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://api.cirrus-ci.com/github/tmigot/RandomKrylov.jl.svg)](https://cirrus-ci.com/github/tmigot/RandomKrylov.jl)
[![Coverage](https://codecov.io/gh/tmigot/RandomKrylov.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/tmigot/RandomKrylov.jl)

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).

## How to Install

Krylov can be installed and tested through the Julia package manager:

```julia
julia> ]
pkg> add RandomKrylov
pkg> test RandomKrylov
```

## Content

This package provides implementations of certain of the most useful Krylov method for a variety of problems:

1. Square or rectangular full-rank systems

<p align="center">
  <b><i>Ax = b</i></b>
</p>

2. Linear least-squares problems

<p align="center">
  minimize ‖<b><i>b</i></b> - <b><i>Ax</i></b>‖
</p>

## References

Some of state-of-art algorithms implemented are presented in the following paper.

Gower, Robert M., and Peter Richtárik. "Randomized iterative methods for linear systems." SIAM Journal on Matrix Analysis and Applications 36.4 (2015): 1660-1690.

We refer to [Krylov.jl](https://github.com/JuliaSmoothOptimizers/Krylov.jl) for (deterministic) Krylov methods in Julia.
This package uses [Stopping.jl](https://github.com/vepiteski/Stopping.jl) as a framework for iterative methods.
