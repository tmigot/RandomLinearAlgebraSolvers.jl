module RandomKrylov

using Krylov

macro kdot(n, x, y)
    return esc(:(Krylov.krylov_dot($n, $x, 1, $y, 1)))
end

using LinearAlgebra, LLSModels, SparseArrays, Stopping
using StatsBase # for sample!

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
include("random_methods.jl")

end
