using RandomLinearAlgebraSolvers
using LinearAlgebra, LLSModels, SparseArrays, Stopping
using Random, Test

using MatrixDepot

Random.seed!(1234)

# algo_list = [:RandomVectorSketch, :randomInv]
algo_list = [:RandomizedKaczmarz, :RandomizedBlockKaczmarz, :RandomizedCD]
algo_pd_list = [:RandomizedCD2, :RandomizedNewton]

@testset "RandomLinearAlgebraSolvers.jl" begin
  for T in [Float16, Float32, Float64], algo in union(algo_list, algo_pd_list)
    A = T[1 0; 0 1]
    b = T[1; 1]
    stp = eval(algo)(RLAStopping(A, b))
    sol = stp.current_state.x
    @test eltype(sol) == T
  end
end

function overdetermined_with_random_matrices(T = Float64, m = 10, n = 5)
  A = rand(T, m, n)
  sol = rand(T, n)
  b = A * sol
  return A, b, sol
end

function overdetermined_with_sparse_random_matrices(T = Float64, m = 10, n = 5)
  A = sprand(T, m, n, 1 / log(n * m))
  sol = rand(T, n)
  b = A * sol
  return A, b, sol
end

function overdetermined_with_random_pd_matrices(T = Float64, n = 10)
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
    stp = RLAStopping(
      A,
      b,
      max_iter = 1000000,
      max_eval = 100000,
      rtol = 1e-2,
      atol = 0.0,
    )
    eval(meth)(stp)
    @test status(stp) == :Optimal
    @show @allocated eval(meth)(stp)
  end
end

for meth in union(algo_list, algo_pd_list), fpb in pb_pd_generators
  @testset "$meth: Test $(fpb)" begin
    A, b, sol = eval(fpb)()
    stp = RLAStopping(
      A,
      b,
      max_iter = 1000000,
      max_eval = 100000,
      rtol = 1e-2,
      atol = 0.0,
    )
    eval(meth)(stp)
    @test status(stp) == :Optimal
    @show @allocated eval(meth)(stp)
  end
end
