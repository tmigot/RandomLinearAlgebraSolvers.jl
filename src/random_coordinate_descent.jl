export RandomizedCD, RandomizedCD2

"""
Randomized coordinate descent

Sect. 3.7 in Gower, R. M., & Richtárik, P. (2015).
Randomized iterative methods for linear systems.
SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.
"""
function RandomizedCD(stp :: AbstractStopping; kwargs...)
    A,b = get_matrix(stp.pb), get_vector(stp.pb)
    m,n = size(A)
    x0  = stp.current_state.x
    xk  = x0
    OK = update_and_start!(stp, x = xk, res = A * xk - b)

    while !OK#norm(res) > prec

     i  = Int(floor(rand() * n)+1) #rand a number between 1 and n
     Ai = A[:,i]

     #ei = zeros(n); ei[i] = 1.0 #unit vector in R^n
     #xk  = Ai == 0 ? x0 : x0 - dot(Ai,res)/norm(Ai,2)^2 * ei
     xk = x0
     xk[i] = x0[i] - dot(Ai, stp.current_state.res)/norm(Ai,2)^2

     OK = update_and_stop!(stp, x = xk, res = A * xk - b)
     x0  = xk

    end

 return stp
end

"""
Randomized coordinate descent
-> for symmetric positive definite matrix

Sect. 3.4 in Gower, R. M., & Richtárik, P. (2015).
Randomized iterative methods for linear systems.
SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.
"""
function RandomizedCD2(stp :: AbstractStopping; kwargs...)
    A,b = get_matrix(stp.pb), get_vector(stp.pb)
    m,n = size(A)
    if (m != n || eigmin(A)<0) throw("RandomizedCD2 error: non-spd matrix") end
    x0  = stp.current_state.x
    xk  = x0
    OK = update_and_start!(stp, x = xk, res = A*xk - b)

    while !OK#norm(res) > prec

     i  = Int(floor(rand() * n)+1) #rand a number between 1 and n
     Ai = A[i,:]

     #ei = zeros(n); ei[i] = 1.0 #unit vector in R^n
     #xk  = Ai == 0 ? x0 : x0 - dot(Ai,res)/norm(Ai,2)^2 * ei
     xk = x0
     xk[i] = A[i,i] == 0 ? x0[i] : x0[i] - (dot(Ai,x0)-b[i])/A[i,i]

     OK = update_and_stop!(stp, x = xk, res = A*xk - b)
     x0  = xk

    end

 return stp
end
