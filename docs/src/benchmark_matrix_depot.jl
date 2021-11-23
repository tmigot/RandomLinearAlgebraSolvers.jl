###############################################################################
#
# Benchmark of methods to solve Ax = b with A an m x n matrix.
# A is generated randomly
# b = A * xref
# x0 = 0
#
# Improvements:
#  * specific symmetric positive definite
#  * use matrices from MatrixDepot
#
###############################################################################
include("methods/linear_algebra_methods.jl")

using MatrixDepot
#Collection of 49 matrices (some are scalable)
#https://github.com/JuliaMatrices/MatrixDepot.jl
#mdinfo()
#listnames(:symmetric)
#mdinfo("hilb")
#A = matrixdepot("hilb", 4)
N = 5
matrix_problems = mdlist(builtin(1:N))
param = 100


#https://juliasmoothoptimizers.github.io/SolverBenchmark.jl/latest/tutorial/
#In this tutorial we illustrate the main uses of SolverBenchmark.
using DataFrames, Printf
#, SolverBenchmark

printstyled("Benchmark linear algebra solvers: \n", color = :green)

#Names of solvers:
names = [:PLSSRRand, :RandomizedBlockKaczmarz]
#problem with: , :RandomVectorSketch,:Juliainv,
#Initialization of the DataFrame for n problems.
stats = Dict(name => DataFrame(
         :id     => 1:N,
         :nvar   => zeros(Int64, N),
         :status => [:Unknown for i = 1:N],
         :time   => NaN*ones(N),
         :iter   => zeros(Int64, N),
         :score  => NaN*ones(N)) for name in names)

for i=1:N

 try
  A = Array(matrixdepot(matrix_problems[i], param)) #convert sparse matrices in normal ones
  xref = 10 * rand(size(A,2))
  b = A * xref

  x0 = zeros(size(A,2))
  la_stop = LinearOperatorStopping(LinearSystem(A, b),
                                   linear_system_check!,
                                   GenericState(x0),
                                   max_iter = 20000,
                                   rtol = sqrt(eps()),
                                   atol = sqrt(eps()))
  @show size(A), cond(Array(A), 2), det(A' * A)
  for name in names

    #solve the problem
    reinit!(la_stop, rstate = true, x = x0)
    la_stop.meta.start_time = time()
    @time eval(name)(la_stop)
    sol = la_stop.current_state.x

    #update the stats from the Stopping
    stats[name].nvar[i]   = size(A,2)
    stats[name].status[i] = status(la_stop)
    stop_has_time = (la_stop.current_state.current_time != nothing)
    stats[name].time[i]   =  stop_has_time ? la_stop.current_state.current_time - la_stop.meta.start_time : time() - la_stop.meta.start_time
    stats[name].iter[i]   = la_stop.meta.nb_of_stop
    stats[name].score[i]  = la_stop.current_state.current_score

  end
catch
end

end #end of main loop

for name in names
    @show name
    @show stats[name]
end

#You can export the table in Latex
#latex_table(stdout, stats[:armijo])

#or run a performance profile:
#using Plots
#pyplot()
#cost(df) = (df.status != :Optimal) * Inf + df.t
#p = performance_profile(stats, cost)
#Plots.svg(p, "profile2")

#or a profile wall:
#solved(df) = (def.status .== :Optimal)
#costs = [df -> .!sovled(df) * Inf + df.t, df -> .!sovled(df) * Inf + df.iter]
#costnames = ["Time", "Iterations"]
#p = profile_solvers(stats, costs, costnames)
#Plots.svg(p, "profile3")

printstyled("The End.", color = :green)