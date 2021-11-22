export RandomizedKaczmarz, RandomizedBlockKaczmarz

"""
Randomized Kaczmarz

Section 3.3 in Gower, R. M., & Richtárik, P. (2015).
Randomized iterative methods for linear systems.
SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.

RK takes a step in the direction of the negative stochastic gradient.
This means that it is equivalent to the SGD method.
However, the stepsize choice is very special: RK chooses the stepsize which leads to the point
which is closest to x* in the Euclidean norm.
"""
function RandomizedKaczmarz(stp::AbstractStopping; is_zero_start::Bool = false, kwargs...)
  A, b = get_matrix(stp.pb), get_vector(stp.pb)
  state = stp.current_state
  m = size(A, 1)
  if is_zero_start
    state.res .= b
  else
    state.res .= b .- A * state.x
  end
  OK = start!(stp, no_opt_check = true)

  while !OK
    i = Int(floor(rand() * m) + 1)
    Ai = view(A, i, :)
    AiAi = dot(Ai, Ai)
    if AiAi != 0
      Aix = dot(Ai, state.x)
      state.x .-= (Aix - b[i]) / AiAi * Ai
    end
    state.res = b - A * state.x
    OK = cheap_stop!(stp)
  end

  return stp
end

"""
Randomized block Kaczmarz

Section 3.5 in Gower, R. M., & Richtárik, P. (2015).
Randomized iterative methods for linear systems.
SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.
"""
function RandomizedBlockKaczmarz(
  stp::AbstractStopping;
  is_zero_start::Bool = false,
  r::Int = 15,
  rand_r::Bool = false,
  kwargs...,
)
  A, b = get_matrix(stp.pb), get_vector(stp.pb)
  state = stp.current_state
  m = size(A, 1)
  if is_zero_start
    state.res .= b
  else
    state.res .= b .- A * state.x
  end
  OK = start!(stp, no_opt_check = true)

  while !OK
    r = rand_r ? Int(floor(rand() * m) + 1) : min(r, m)
    sub = zeros(Int64, r)
    #sample!(1:m, sub) #x is a subset of [1,...,m] of size r #sample repeat entries
    sub = StatsBase.randperm(m)[1:r]

    Ai = view(A, sub, :)
    resk = state.res[sub]
    if !(Ai == 0)
      state.x .-= Ai' * pinv(Matrix(Ai * Ai')) * resk
    end
    state.res .= A * state.x .- b
    OK = cheap_stop!(stp)
  end

  return stp
end
