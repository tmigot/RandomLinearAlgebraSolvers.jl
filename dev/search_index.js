var documenterSearchIndex = {"docs":
[{"location":"matrixdepot/#Benchmark-on-MatrixDepot-matrices","page":"Benchmark MatrixDepot","title":"Benchmark on MatrixDepot matrices","text":"","category":"section"},{"location":"matrixdepot/","page":"Benchmark MatrixDepot","title":"Benchmark MatrixDepot","text":"using RandomLinearAlgebraSolvers, Random, LinearAlgebra, SparseArrays\nRandom.seed!(1234)","category":"page"},{"location":"matrixdepot/","page":"Benchmark MatrixDepot","title":"Benchmark MatrixDepot","text":"We use here MatrixDepot a collection of 49 matrices (some are scalable).","category":"page"},{"location":"matrixdepot/","page":"Benchmark MatrixDepot","title":"Benchmark MatrixDepot","text":"using MatrixDepot","category":"page"},{"location":"matrixdepot/","page":"Benchmark MatrixDepot","title":"Benchmark MatrixDepot","text":"Benchmark of methods to solve Ax = b with  b = A * xref.","category":"page"},{"location":"matrixdepot/","page":"Benchmark MatrixDepot","title":"Benchmark MatrixDepot","text":"using DataFrames, Printf, SolverBenchmark, Stopping","category":"page"},{"location":"matrixdepot/","page":"Benchmark MatrixDepot","title":"Benchmark MatrixDepot","text":"N = 40 # number of problems\nmatrix_problems = mdlist(builtin(1:N))\nparam = 50","category":"page"},{"location":"matrixdepot/","page":"Benchmark MatrixDepot","title":"Benchmark MatrixDepot","text":"Names of solvers:","category":"page"},{"location":"matrixdepot/","page":"Benchmark MatrixDepot","title":"Benchmark MatrixDepot","text":"names = [:RandomizedKaczmarz, :RandomizedBlockKaczmarz, :RandomizedCD]","category":"page"},{"location":"matrixdepot/","page":"Benchmark MatrixDepot","title":"Benchmark MatrixDepot","text":"#Initialization of the DataFrame for n problems.\nstats = Dict(name => DataFrame(\n         :id     => 1:N,\n         :nvar   => zeros(Int64, N),\n         :status => [:Unknown for i = 1:N],\n         :time   => NaN*ones(N),\n         :iter   => zeros(Int64, N),\n         :score  => NaN*ones(N)) for name in names)","category":"page"},{"location":"matrixdepot/","page":"Benchmark MatrixDepot","title":"Benchmark MatrixDepot","text":"for i=1:N\n\n  A = matrixdepot(matrix_problems[i], param)\n  m, n = size(A)\n  xref = 100 * rand(n)\n  b = A * xref\n\n  x0 = zeros(size(A,2))\n  la_stop = LAStopping(A, b, max_iter = 100000, rtol = sqrt(eps()), atol = sqrt(eps()))\n  for name in names\n\n    #solve the problem\n    reinit!(la_stop, rstate = true, x = x0, res = similar(b))\n    la_stop.meta.start_time = time()\n    @time eval(name)(la_stop, r = 80, is_zero_start=true)\n    sol = la_stop.current_state.x\n\n    #update the stats from the Stopping\n    stats[name].nvar[i]   = n\n    stats[name].status[i] = status(la_stop)\n    stop_has_time = (la_stop.current_state.current_time != nothing)\n    stats[name].time[i]   =  stop_has_time ? la_stop.current_state.current_time - la_stop.meta.start_time : time() - la_stop.meta.start_time\n    stats[name].iter[i]   = la_stop.meta.nb_of_stop\n    stats[name].score[i]  = norm(la_stop.current_state.current_score, Inf)\n\n  end\n\nend","category":"page"},{"location":"matrixdepot/","page":"Benchmark MatrixDepot","title":"Benchmark MatrixDepot","text":"for name in names\n    @show name\n    @show stats[name]\nend","category":"page"},{"location":"matrixdepot/","page":"Benchmark MatrixDepot","title":"Benchmark MatrixDepot","text":"or run a performance profile:","category":"page"},{"location":"matrixdepot/","page":"Benchmark MatrixDepot","title":"Benchmark MatrixDepot","text":"using Plots\ngr()\ncost(df) = (df.status .!= :Optimal) * Inf + df.time\np = performance_profile(stats, cost)","category":"page"},{"location":"sdp/#Benchmark-on-dense-random-symmetric-positive-definite-matrices","page":"Benchmark sdp matrices","title":"Benchmark on dense random symmetric positive definite matrices","text":"","category":"section"},{"location":"sdp/","page":"Benchmark sdp matrices","title":"Benchmark sdp matrices","text":"using RandomLinearAlgebraSolvers, Random, LinearAlgebra, SparseArrays\nRandom.seed!(1234)","category":"page"},{"location":"sdp/","page":"Benchmark sdp matrices","title":"Benchmark sdp matrices","text":"Benchmark of methods to solve Ax = b with A an n x n randomly generated symmetric positive definite matrix, and b = A * xref.","category":"page"},{"location":"sdp/","page":"Benchmark sdp matrices","title":"Benchmark sdp matrices","text":"using DataFrames, Printf, SolverBenchmark, Stopping","category":"page"},{"location":"sdp/","page":"Benchmark sdp matrices","title":"Benchmark sdp matrices","text":"N = 5 # number of problems\nn = 100 # size of A: n x n","category":"page"},{"location":"sdp/","page":"Benchmark sdp matrices","title":"Benchmark sdp matrices","text":"Names of solvers:","category":"page"},{"location":"sdp/","page":"Benchmark sdp matrices","title":"Benchmark sdp matrices","text":"names = [:RandomizedKaczmarz, :RandomizedBlockKaczmarz, :RandomizedCD, :RandomizedCD2, :RandomizedNewton]","category":"page"},{"location":"sdp/","page":"Benchmark sdp matrices","title":"Benchmark sdp matrices","text":"#Initialization of the DataFrame for n problems.\nstats = Dict(name => DataFrame(\n         :id     => 1:N,\n         :nvar   => zeros(Int64, N),\n         :status => [:Unknown for i = 1:N],\n         :time   => NaN*ones(N),\n         :iter   => zeros(Int64, N),\n         :score  => NaN*ones(N)) for name in names)","category":"page"},{"location":"sdp/","page":"Benchmark sdp matrices","title":"Benchmark sdp matrices","text":"for i=1:N\n\n  Ar = rand(n,n); A = 1/sqrt(n) * Ar'*Ar + Matrix{Float64}(I, n, n)\n  xref = 100 * rand(n)\n  b = A * xref\n\n  x0 = zeros(size(A,2))\n  la_stop = LAStopping(A, b, max_iter = 100000, rtol = sqrt(eps()), atol = sqrt(eps()))\n  for name in names\n\n    #solve the problem\n    reinit!(la_stop, rstate = true, x = x0, res = similar(b))\n    la_stop.meta.start_time = time()\n    @time eval(name)(la_stop, r = 80, is_zero_start=true)\n    sol = la_stop.current_state.x\n\n    #update the stats from the Stopping\n    stats[name].nvar[i]   = n\n    stats[name].status[i] = status(la_stop)\n    stop_has_time = (la_stop.current_state.current_time != nothing)\n    stats[name].time[i]   =  stop_has_time ? la_stop.current_state.current_time - la_stop.meta.start_time : time() - la_stop.meta.start_time\n    stats[name].iter[i]   = la_stop.meta.nb_of_stop\n    stats[name].score[i]  = norm(la_stop.current_state.current_score, Inf)\n\n  end\n\nend","category":"page"},{"location":"sdp/","page":"Benchmark sdp matrices","title":"Benchmark sdp matrices","text":"for name in names\n    @show name\n    @show stats[name]\nend","category":"page"},{"location":"sdp/","page":"Benchmark sdp matrices","title":"Benchmark sdp matrices","text":"or run a performance profile:","category":"page"},{"location":"sdp/","page":"Benchmark sdp matrices","title":"Benchmark sdp matrices","text":"using Plots\ngr()\ncost(df) = (df.status .!= :Optimal) * Inf + df.time\np = performance_profile(stats, cost)","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = RandomLinearAlgebraSolvers","category":"page"},{"location":"#RandomLinearAlgebraSolvers","page":"Home","title":"RandomLinearAlgebraSolvers","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for RandomLinearAlgebraSolvers.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [RandomLinearAlgebraSolvers]","category":"page"},{"location":"#RandomLinearAlgebraSolvers.RLAStopping-Union{Tuple{S}, Tuple{Any, S}} where S","page":"Home","title":"RandomLinearAlgebraSolvers.RLAStopping","text":"stp = RLAStopping(A, b::S; n_listofstates::Int = 0, kwargs...)\n\nCreator of a LAStopping tailored for this package. The problem, stp.pb, is a Stopping.LinearSystem if A is dense, and a LLSModels.LLSModel otherwise. The state allocates space for the residual, stp.current_state.res, of length |b|. This stopping uses the infinity norm of stp.current_state.res to declare optimality.\n\n\n\n\n\n","category":"method"},{"location":"#RandomLinearAlgebraSolvers.RandomVectorSketch-Tuple{Stopping.AbstractStopping}","page":"Home","title":"RandomLinearAlgebraSolvers.RandomVectorSketch","text":"Random vector sketch\n\nSection 3.2 in Gower, R. M., & Richtárik, P. (2015). Randomized iterative methods for linear systems. SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.\n\n\n\n\n\n","category":"method"},{"location":"#RandomLinearAlgebraSolvers.RandomizedBlockKaczmarz-Tuple{Stopping.AbstractStopping}","page":"Home","title":"RandomLinearAlgebraSolvers.RandomizedBlockKaczmarz","text":"Randomized block Kaczmarz\n\nSection 3.5 in Gower, R. M., & Richtárik, P. (2015). Randomized iterative methods for linear systems. SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.\n\n\n\n\n\n","category":"method"},{"location":"#RandomLinearAlgebraSolvers.RandomizedCD-Tuple{Stopping.AbstractStopping}","page":"Home","title":"RandomLinearAlgebraSolvers.RandomizedCD","text":"Randomized coordinate descent\n\nSection 3.7 in Gower, R. M., & Richtárik, P. (2015). Randomized iterative methods for linear systems. SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.\n\n\n\n\n\n","category":"method"},{"location":"#RandomLinearAlgebraSolvers.RandomizedCD2-Tuple{Stopping.AbstractStopping}","page":"Home","title":"RandomLinearAlgebraSolvers.RandomizedCD2","text":"Randomized coordinate descent for symmetric positive definite matrix\n\nSection 3.4 in Gower, R. M., & Richtárik, P. (2015). Randomized iterative methods for linear systems. SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.\n\n\n\n\n\n","category":"method"},{"location":"#RandomLinearAlgebraSolvers.RandomizedKaczmarz-Tuple{Stopping.AbstractStopping}","page":"Home","title":"RandomLinearAlgebraSolvers.RandomizedKaczmarz","text":"Randomized Kaczmarz\n\nSection 3.3 in Gower, R. M., & Richtárik, P. (2015). Randomized iterative methods for linear systems. SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.\n\nRK takes a step in the direction of the negative stochastic gradient. This means that it is equivalent to the SGD method. However, the stepsize choice is very special: RK chooses the stepsize which leads to the point which is closest to x* in the Euclidean norm.\n\n\n\n\n\n","category":"method"},{"location":"#RandomLinearAlgebraSolvers.RandomizedNewton-Tuple{Stopping.AbstractStopping}","page":"Home","title":"RandomLinearAlgebraSolvers.RandomizedNewton","text":"Randomized Newton -> for symmetric positive definite matrix\n\nSection 3.6 in Gower, R. M., & Richtárik, P. (2015). Randomized iterative methods for linear systems. SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.\n\n\n\n\n\n","category":"method"},{"location":"#RandomLinearAlgebraSolvers.randomInv-Tuple{Stopping.AbstractStopping}","page":"Home","title":"RandomLinearAlgebraSolvers.randomInv","text":"Solve TAx=Tb with T a random projector\n\n\n\n\n\n","category":"method"},{"location":"#RandomLinearAlgebraSolvers.random_matrix_1-Tuple{Any, Any}","page":"Home","title":"RandomLinearAlgebraSolvers.random_matrix_1","text":"Random projector of size k x m with a normal distribution\n\n\n\n\n\n","category":"method"},{"location":"#RandomLinearAlgebraSolvers.random_matrix_2-Tuple{Any, Any}","page":"Home","title":"RandomLinearAlgebraSolvers.random_matrix_2","text":"Random projector of size k x m with -1 or 1 both with probability 1/2\n\n\n\n\n\n","category":"method"},{"location":"#RandomLinearAlgebraSolvers.random_matrix_3-Tuple{Any, Any}","page":"Home","title":"RandomLinearAlgebraSolvers.random_matrix_3","text":"Random projector of size k x m with -1,0,1 respectively with probability 1/6,4/6,1/6\n\n\n\n\n\n","category":"method"},{"location":"#RandomLinearAlgebraSolvers.random_matrix_4-Tuple{Any, Any}","page":"Home","title":"RandomLinearAlgebraSolvers.random_matrix_4","text":"Random projector of size k x m with orthogonal projection on a random k-dimensional linear subspace of R^m\n\n\n\n\n\n","category":"method"},{"location":"#RandomLinearAlgebraSolvers.random_projector_rank-Tuple{Function, Int64, Int64, Int64}","page":"Home","title":"RandomLinearAlgebraSolvers.random_projector_rank","text":"Check the average rank of random projector over N random matrices of size k*m\n\n\n\n\n\n","category":"method"},{"location":"#RandomLinearAlgebraSolvers.random_projector_sparse-Tuple{Function, Int64, Int64, Int64}","page":"Home","title":"RandomLinearAlgebraSolvers.random_projector_sparse","text":"Check the average sparsity of random projector over N random matrices of size k*m\n\n\n\n\n\n","category":"method"},{"location":"rectangular/#Benchmark-on-dense-overdetermined-random-matrices","page":"Benchmark overdetermined","title":"Benchmark on dense overdetermined random matrices","text":"","category":"section"},{"location":"rectangular/","page":"Benchmark overdetermined","title":"Benchmark overdetermined","text":"using RandomLinearAlgebraSolvers, Random, LinearAlgebra, SparseArrays\nRandom.seed!(1234)","category":"page"},{"location":"rectangular/","page":"Benchmark overdetermined","title":"Benchmark overdetermined","text":"Benchmark of methods to solve Ax = b with A a randomly generated m x n matrix, and b = A * xref.","category":"page"},{"location":"rectangular/","page":"Benchmark overdetermined","title":"Benchmark overdetermined","text":"using DataFrames, Printf, SolverBenchmark, Stopping","category":"page"},{"location":"rectangular/","page":"Benchmark overdetermined","title":"Benchmark overdetermined","text":"N = 5 # number of problems\nm, n = 10000, 100 # size of A: m x n","category":"page"},{"location":"rectangular/","page":"Benchmark overdetermined","title":"Benchmark overdetermined","text":"Names of solvers:","category":"page"},{"location":"rectangular/","page":"Benchmark overdetermined","title":"Benchmark overdetermined","text":"names = [:RandomizedKaczmarz, :RandomizedBlockKaczmarz, :RandomizedCD]","category":"page"},{"location":"rectangular/","page":"Benchmark overdetermined","title":"Benchmark overdetermined","text":"#Initialization of the DataFrame for n problems.\nstats = Dict(name => DataFrame(\n         :id     => 1:N,\n         :nvar   => zeros(Int64, N),\n         :status => [:Unknown for i = 1:N],\n         :time   => NaN*ones(N),\n         :iter   => zeros(Int64, N),\n         :score  => NaN*ones(N)) for name in names)","category":"page"},{"location":"rectangular/","page":"Benchmark overdetermined","title":"Benchmark overdetermined","text":"for i=1:N\n\n  A = 100 * rand(m, n)\n  xref = 100 * rand(n)\n  b = A * xref\n\n  x0 = zeros(size(A,2))\n  la_stop = LAStopping(A, b, max_iter = 100000, rtol = sqrt(eps()), atol = sqrt(eps()))\n  for name in names\n\n    #solve the problem\n    reinit!(la_stop, rstate = true, x = x0, res = similar(b))\n    la_stop.meta.start_time = time()\n    @time eval(name)(la_stop, r = 80, is_zero_start=true)\n    sol = la_stop.current_state.x\n\n    #update the stats from the Stopping\n    stats[name].nvar[i]   = n\n    stats[name].status[i] = status(la_stop)\n    stop_has_time = (la_stop.current_state.current_time != nothing)\n    stats[name].time[i]   =  stop_has_time ? la_stop.current_state.current_time - la_stop.meta.start_time : time() - la_stop.meta.start_time\n    stats[name].iter[i]   = la_stop.meta.nb_of_stop\n    stats[name].score[i]  = norm(la_stop.current_state.current_score, Inf)\n\n  end\n\nend","category":"page"},{"location":"rectangular/","page":"Benchmark overdetermined","title":"Benchmark overdetermined","text":"for name in names\n    @show name\n    @show stats[name]\nend","category":"page"},{"location":"rectangular/","page":"Benchmark overdetermined","title":"Benchmark overdetermined","text":"or run a performance profile:","category":"page"},{"location":"rectangular/","page":"Benchmark overdetermined","title":"Benchmark overdetermined","text":"using Plots\ngr()\ncost(df) = (df.status .!= :Optimal) * Inf + df.time\np = performance_profile(stats, cost)","category":"page"}]
}
