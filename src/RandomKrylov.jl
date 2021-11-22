module RandomKrylov

using LinearAlgebra, LLSModels, SparseArrays, Stopping
using StatsBase # for sample!

@inline krylov_dot(n :: Integer, x :: Vector{T}, dx :: Integer, y :: Vector{T}, dy :: Integer) where T <: BLAS.BlasReal = BLAS.dot(n, x, dx, y, dy)
@inline krylov_dot(n :: Integer, x :: Vector{T}, dx :: Integer, y :: Vector{T}, dy :: Integer) where T <: BLAS.BlasComplex = BLAS.dotc(n, x, dx, y, dy)
@inline krylov_dot(n :: Integer, x :: AbstractVector{T}, dx :: Integer, y :: AbstractVector{T}, dy :: Integer) where T <: Number = dot(x, y)

macro kdot(n, x, y)
  return esc(:(krylov_dot($n, $x, 1, $y, 1)))
end

function get_matrix(lls::LLSModel)
  return sparse(lls.Arows, lls.Acols, lls.Avals, length(lls.b), lls.meta.nvar)
end

function get_vector(lls::LLSModel)
  return lls.b
end

function get_matrix(pb::LinearSystem)
  return pb.A
end

function get_vector(pb::LinearSystem)
  return pb.b
end

include("random_Kaczmarcz.jl")
include("random_coordinate_descent.jl")
include("random_projector.jl")
include("random_Newton.jl")
include("random_methods.jl")

end
