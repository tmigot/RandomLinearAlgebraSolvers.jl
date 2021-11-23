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

This package provides implementations of certain of the randomized numerical methods for linear algebra:
<p align="center">
  minimize ‖<b><i>b</i></b> - <b><i>Ax</i></b>‖
</p>
It includes classical Kaczmarcz and coordinate descent methods.

The package also provides random projectors used to solve
<p align="center">
  minimize ‖<b><i>Tb</i></b> - <b><i>TAx</i></b>‖
</p>
as an approximation of the initial system in the sense of the [Johnson–Lindenstrauss lemma](https://en.wikipedia.org/wiki/Johnson–Lindenstrauss_lemma).

## Example

This package uses [Stopping.jl](https://github.com/vepiteski/Stopping.jl) as a framework for iterative methods.
```
using RandomKrylov
A, b = rand(10, 5), rand(10)
stp = RLAStopping(A, b)
RandomizedKaczmarz(stp)
```

## References

Some of state-of-art algorithms implemented are presented in the following paper.

> Gower, Robert M., and Peter Richtárik. "Randomized iterative methods for linear systems." SIAM Journal on Matrix Analysis and Applications 36.4 (2015): 1660-1690.

We refer to [Krylov.jl](https://github.com/JuliaSmoothOptimizers/Krylov.jl) for deterministic Krylov methods in Julia.

# Bug reports and discussions

If you think you found a bug, feel free to open an [issue](https://github.com/tmigot/RandomKrylov.jl/issues).
Focused suggestions and requests can also be opened as issues. Before opening a pull request, start an issue or a discussion on the topic, please.

The package is still at an early stage and new contributions are very welcome. We would like to gather as much information as possible on the provenance of new algorithms.
