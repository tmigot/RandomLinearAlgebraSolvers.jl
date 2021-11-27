## Benchmark on dense random symmetric positive definite matrices

```@example ex1
using RandomLinearAlgebraSolvers, Random, LinearAlgebra, SparseArrays
Random.seed!(1234)
```
Benchmark of methods to solve `Ax = b` with `A` an `n x n` randomly generated symmetric positive definite matrix, and `b = A * xref`.
```@example ex1
using DataFrames, Printf, SolverBenchmark, Stopping
```
```@example ex1
N = 5 # number of problems
n = 100 # size of A: n x n
```
Names of solvers:
```@example ex1
names = [:RandomizedKaczmarz, :RandomizedBlockKaczmarz, :RandomizedCD, :RandomizedCD2, :RandomizedNewton]
```

```@example ex1
#Initialization of the DataFrame for n problems.
stats = Dict(name => DataFrame(
         :id     => 1:N,
         :nvar   => zeros(Int64, N),
         :status => [:Unknown for i = 1:N],
         :time   => NaN*ones(N),
         :iter   => zeros(Int64, N),
         :score  => NaN*ones(N)) for name in names)
```

```@example ex1
for i=1:N

  Ar = rand(n,n); A = 1/sqrt(n) * Ar'*Ar + Matrix{Float64}(I, n, n)
  xref = 100 * rand(n)
  b = A * xref

  x0 = zeros(size(A,2))
  la_stop = RLAStopping(A, b, max_iter = 100000, rtol = sqrt(eps()), atol = sqrt(eps()))
  for name in names

    #solve the problem
    reinit!(la_stop, rstate = true, x = x0, res = similar(b))
    la_stop.meta.start_time = time()
    @time eval(name)(la_stop, r = 80, is_zero_start=true)
    sol = la_stop.current_state.x

    #update the stats from the Stopping
    stats[name].nvar[i]   = n
    stats[name].status[i] = status(la_stop)
    stop_has_time = (la_stop.current_state.current_time != nothing)
    stats[name].time[i]   =  stop_has_time ? la_stop.current_state.current_time - la_stop.meta.start_time : time() - la_stop.meta.start_time
    stats[name].iter[i]   = la_stop.meta.nb_of_stop
    stats[name].score[i]  = norm(la_stop.current_state.current_score, Inf)

  end

end
```

```@example ex1
for name in names
    @show name
    @show stats[name]
end
```
or run a performance profile:
```@example ex1
using Plots
gr()
cost(df) = (df.status .!= :Optimal) * Inf + df.time
p = performance_profile(stats, cost)
```
