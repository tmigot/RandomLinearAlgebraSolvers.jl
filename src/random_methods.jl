export randomInv, RandomVectorSketch, RandomizedNewton

"""
Solve TAx=Tb with T a random projector
"""
function randomInv(stp :: AbstractStopping; is_zero_start::Bool = false, k = 50, random_func = random_matrix_1, kwargs...)
  A, b = get_matrix(stp.pb), get_vector(stp.pb)
  m = size(A, 1)

  start!(stp)
  k   = Int(floor(m/8))
  T   = random_func(k, m)
  xk =  T * A \ T * b

  return update_and_stop!(stp, x = xk)
end

"""
Random vector sketch

Sect. 3.2 in Gower, R. M., & Richtárik, P. (2015).
Randomized iterative methods for linear systems.
SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.
"""
function RandomVectorSketch(stp :: AbstractStopping;
    is_zero_start::Bool = false, random_func = random_matrix_1, kwargs...)
 A,b = get_matrix(stp.pb), get_vector(stp.pb)
 m,n = size(A)
 x0  = stp.current_state.x
 xk  = x0
 stp.current_state.res = A*xk - b
 OK = start!(stp)

 while !OK #norm(res) > prec

  s  = vec(random_func(1, m)')
  As  = A'*s
  xk  = As == 0 ? x0 : x0 - dot(s,stp.current_state.res)/dot(As,As) * As

  update!(stp.current_state, x = xk)
  stp.current_state.res = A*xk - b
  OK = stop!(stp)
  x0  = xk

 end

 return stp
end

"""
Randomized Newton
-> for symmetric positive definite matrix

Sect. 3.6 in Gower, R. M., & Richtárik, P. (2015).
Randomized iterative methods for linear systems.
SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.
"""
function RandomizedNewton(stp :: AbstractStopping; is_zero_start::Bool = false, r :: Int = 15, rand_r :: Bool = false, kwargs...)
    A,b = get_matrix(stp.pb), get_vector(stp.pb)
    m,n = size(A)
    if (m != n || eigmin(A)<0) throw("RandomizedNewton error: non-spd matrix") end
    x0  = stp.current_state.x
    xk  = x0
    stp.current_state.res = A*xk - b
    OK = start!(stp)

    while !OK#norm(res) > prec

     #r   = Int(floor(rand() * n)+1) #rand a number between 1 and n
     r    = rand_r ? Int(floor(rand() * n)+1) : r
     sub = zeros(Int64, r)
     sample!(1:n, sub) #x is a subset of [1,...,n] of size r

     Ai = A[sub,:]
     resk = stp.current_state.res[sub] #resk = Ai*xk - b[sub]

     try
      d  = Ai \ -resk
      xk  = x0 + d
    catch #if Ai d = -resk has no solution
      xk = x0
     end

     update!(stp.current_state, x = xk)
     stp.current_state.res = A*xk - b
     OK = stop!(stp)
     x0  = xk

    end

 return stp
end
