export randomInv, RandomVectorSketch

"""
Solve TAx=Tb with T a random projector
"""
function randomInv(
  stp::AbstractStopping;
  is_zero_start::Bool = false,
  k = 50,
  random_func = random_matrix_1,
  kwargs...,
)
  A, b = get_matrix(stp.pb), get_vector(stp.pb)
  m = size(A, 1)

  start!(stp)
  k = Int(floor(m / 8))
  T = random_func(k, m)
  xk = T * A \ T * b

  return update_and_stop!(stp, x = xk)
end

"""
Random vector sketch

Sect. 3.2 in Gower, R. M., & Richtárik, P. (2015).
Randomized iterative methods for linear systems.
SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.
"""
function RandomVectorSketch(
  stp::AbstractStopping;
  is_zero_start::Bool = false,
  random_func = random_matrix_1,
  kwargs...,
)
  A, b = get_matrix(stp.pb), get_vector(stp.pb)
  m, n = size(A)
  x0 = stp.current_state.x
  xk = x0
  stp.current_state.res = A * xk - b
  OK = start!(stp)

  while !OK
    s = vec(random_func(1, m)')
    As = A' * s
    xk = As == 0 ? x0 : x0 - dot(s, stp.current_state.res) / dot(As, As) * As

    update!(stp.current_state, x = xk)
    stp.current_state.res = A * xk - b
    OK = stop!(stp)
    x0 = xk
  end

  return stp
end
