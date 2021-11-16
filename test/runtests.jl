using RandomKrylov
using LinearAlgebra, LLSModels, SparseArrays, Stopping
using Random, Test

using MatrixDepot

Random.seed!(1234)

# algo_list = [:RandomVectorSketch, :randomInv]
algo_list = [:RandomizedKaczmarz, :RandomizedBlockKaczmarz, :RandomizedCD]
algo_pd_list = [:RandomizedCD2, :RandomizedNewton]

@testset "RandomKrylov.jl" begin
  # Write your tests here.
end

function overdetermined_with_random_matrices(T = Float64, m = 10, n = 5) # , m=1000, n=500
  A = rand(T, m, n)
  sol = rand(T, n)
  b = A * sol
  return A, b, sol
end

function overdetermined_with_sparse_random_matrices(T = Float64, m = 10, n = 5) # m=1000, n=500
  A = sprand(T, m, n, 1 / log(n * m))
  sol = rand(T, n)
  b = A * sol
  return A, b, sol
end

function overdetermined_with_random_pd_matrices(T = Float64, n = 10) # n=1000
  A = rand(T, n, n)
  A = A' * A
  A = (A + A') / 2
  sol = rand(T, n)
  b = A * sol
  return A, b, sol
end

function hilbert(T = Float64, n = 4)
  A = matrixdepot("hilb", n)
  sol = rand(T, n)
  b = A * sol
  return A, b, sol
end

pb_generators =
  [:hilbert, :overdetermined_with_random_matrices, :overdetermined_with_sparse_random_matrices]

pb_pd_generators = [:overdetermined_with_random_pd_matrices]

for meth in algo_list, fpb in pb_generators
  @testset "$meth: Test $(fpb)" begin
    A, b, sol = eval(fpb)()
    pb = issparse(A) ? LLSModel(A, b) : LinearSystem(A, b)
    state = GenericState(similar(sol), similar(b), res = similar(b))
    stp = LAStopping(
      pb,
      state,
      max_iter = 1000000,
      max_eval = 100000,
      rtol = 1e-2,
      atol = 0.0,
      optimality_check = (pb, state) -> state.res,
    )
    eval(meth)(stp)
    @test status(stp) == :Optimal
    @show @allocated eval(meth)(stp)
  end
end

for meth in union(algo_list, algo_pd_list), fpb in pb_pd_generators
  @testset "$meth: Test $(fpb)" begin
    A, b, sol = eval(fpb)()
    pb = issparse(A) ? LLSModel(A, b) : LinearSystem(A, b)
    state = GenericState(similar(sol), similar(b), res = similar(b))
    stp = LAStopping(
      pb,
      state,
      max_iter = 1000000,
      max_eval = 100000,
      rtol = 1e-2,
      atol = 0.0,
      optimality_check = (pb, state) -> state.res,
    )
    eval(meth)(stp)
    @test status(stp) == :Optimal
    @show @allocated eval(meth)(stp)
  end
end
