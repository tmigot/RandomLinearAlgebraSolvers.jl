export RandomizedNewton

"""
Randomized Newton
-> for symmetric positive definite matrix

Sect. 3.6 in Gower, R. M., & Richtárik, P. (2015).
Randomized iterative methods for linear systems.
SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.
"""
function RandomizedNewton(
  stp::AbstractStopping;
  is_zero_start::Bool = false,
  r::Int = 15,
  rand_r::Bool = false,
  kwargs...,
)
  A, b = get_matrix(stp.pb), get_vector(stp.pb)
  m, n = size(A)
  state = stp.current_state
  if is_zero_start
    state.res .= b
  else
    state.res .= b .- A * state.x
  end
  OK = start!(stp, no_opt_check = true)

  while !OK

    #r   = Int(floor(rand() * n)+1) #rand a number between 1 and n
    r = rand_r ? Int(floor(rand() * n) + 1) : r
    sub = zeros(Int64, r)
    sample!(1:n, sub) #x is a subset of [1,...,n] of size r

    Ai = A[sub, :]
    resk = state.res[sub] #resk = Ai*xk - b[sub]

    if 0 ∉ resk
      state.x .-= Ai \ resk
    end

    state.res .= b .- A * state.x
    OK = cheap_stop!(stp)
  end

  return stp
end
