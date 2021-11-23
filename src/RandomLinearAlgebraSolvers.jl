module RandomLinearAlgebraSolvers

using LinearAlgebra, LLSModels, SparseArrays, Stopping
using StatsBase # for sample!

@inline krylov_dot(
  n::Integer,
  x::Vector{T},
  dx::Integer,
  y::Vector{T},
  dy::Integer,
) where {T <: BLAS.BlasReal} = BLAS.dot(n, x, dx, y, dy)
@inline krylov_dot(
  n::Integer,
  x::Vector{T},
  dx::Integer,
  y::Vector{T},
  dy::Integer,
) where {T <: BLAS.BlasComplex} = BLAS.dotc(n, x, dx, y, dy)
@inline krylov_dot(
  n::Integer,
  x::AbstractVector{T},
  dx::Integer,
  y::AbstractVector{T},
  dy::Integer,
) where {T <: Number} = dot(x, y)

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

export RLAStopping

"""
    stp = RLAStopping(A, b::S; n_listofstates::Int = 0, kwargs...)

Creator of a `LAStopping` tailored for this package.
The problem, stp.pb, is a `Stopping.LinearSystem` if `A` is dense, and a `LLSModels.LLSModel` otherwise.
The state allocates space for the residual, `stp.current_state.res`, of length `|b|`.
This stopping uses the infinity norm of `stp.current_state.res` to declare optimality.
"""
function RLAStopping(A, b::S; n_listofstates::Int = 0, kwargs...) where {S}
  T = eltype(S)
  sparse = issparse(A)
  pb = sparse ? LLSModel(A, b) : LinearSystem(A, b)
  x0 = zeros(T, size(A, 2))
  state = GenericState(similar(x0), similar(b), res = similar(b))

  mcntrs = sparse ? init_max_counters_NLS() : init_max_counters_linear_operators()

  if n_listofstates > 0 && :list ∉ keys(kwargs)
    list = ListofStates(n_listofstates, Val{typeof(state)}())
    return LAStopping(
      pb,
      state,
      max_cntrs = mcntrs,
      list = list,
      optimality_check = (pb, state) -> state.res;
      kwargs...,
    )
  end

  return LAStopping(
    pb,
    state,
    max_cntrs = mcntrs,
    optimality_check = (pb, state) -> state.res;
    kwargs...,
  )
end

include("random_Kaczmarcz.jl")
include("random_coordinate_descent.jl")
include("random_projector.jl")
include("random_Newton.jl")
include("random_methods.jl")

end
