export RandomizedCD, RandomizedCD2

"""
Randomized coordinate descent

Section 3.7 in Gower, R. M., & Richtárik, P. (2015).
Randomized iterative methods for linear systems.
SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.
"""
function RandomizedCD(stp::AbstractStopping; is_zero_start::Bool = false, kwargs...)
  A, b = get_matrix(stp.pb), get_vector(stp.pb)
  state = stp.current_state
  m, n = size(A)
  if is_zero_start
    state.res .= b
  else
    state.res .= b .- A * state.x
  end
  OK = start!(stp, no_opt_check = true)

  while !OK
    i = Int(floor(rand() * n) + 1)
    Ai = view(A, :, i)

    Aires = @kdot(m, Ai, state.res)
    nAi = @kdot(m, Ai, Ai)
    state.x[i] += Aires / nAi
    state.res .-= Ai .* (Aires / nAi)
    OK = cheap_stop!(stp)
  end

  return stp
end

"""
Randomized coordinate descent for symmetric positive definite matrix

Section 3.4 in Gower, R. M., & Richtárik, P. (2015).
Randomized iterative methods for linear systems.
SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.
"""
function RandomizedCD2(stp::AbstractStopping; is_zero_start::Bool = false, kwargs...)
  A, b = get_matrix(stp.pb), get_vector(stp.pb)
  state = stp.current_state
  n = size(A, 2)
  if is_zero_start
    state.res .= b
  else
    state.res .= b .- A * state.x
  end
  OK = start!(stp, no_opt_check = true)

  while !OK
    i = Int(floor(rand() * n) + 1)
    Ai = view(A, i, :)

    Aix = @kdot(n, Ai, state.x)
    if A[i, i] != 0
      state.x[i] -= (Aix - b[i]) / A[i, i]
    end
    state.res .+= Ai .* ((Aix - b[i]) / A[i, i])

    OK = cheap_stop!(stp)
  end

  return stp
end
