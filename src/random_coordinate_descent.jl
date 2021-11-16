export RandomizedCD, RandomizedCD2

"""
Randomized coordinate descent

Sect. 3.7 in Gower, R. M., & Richtárik, P. (2015).
Randomized iterative methods for linear systems.
SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.
"""
function RandomizedCD(
    stp :: AbstractStopping;
    is_zero_start::Bool = false,
    kwargs...)
    A, b = get_matrix(stp.pb), get_vector(stp.pb)
    state = stp.current_state
    m, n = size(A)
    x  = state.x
    state.res = is_zero_start ? b : b - A * x
    OK = start!(stp, no_start_opt_check = true)

    while !OK

     i = Int(floor(rand() * n)+1) #rand a number between 1 and n
     Ai = view(A, :,i)

     #ei = zeros(n); ei[i] = 1.0 #unit vector in R^n
     #xk  = Ai == 0 ? x0 : x0 - dot(Ai,res)/norm(Ai,2)^2 * ei
     Aires = @kdot(m, Ai, state.res)
     nAi = @kdot(m, Ai, Ai)
     state.x[i] += Aires/nAi
     state.res .-= Ai .* (Aires / nAi)
     OK = cheap_stop!(stp)
    end

 return stp
end

"""
Randomized coordinate descent for symmetric positive definite matrix

Sect. 3.4 in Gower, R. M., & Richtárik, P. (2015).
Randomized iterative methods for linear systems.
SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.
"""
function RandomizedCD2(stp :: AbstractStopping; is_zero_start::Bool = false, kwargs...)
    A,b = get_matrix(stp.pb), get_vector(stp.pb)
    state = stp.current_state
    n = size(A, 2)
    x  = state.x
    state.res = is_zero_start ? b : b - A * x
    OK = start!(stp, no_start_opt_check = true)

    #if (m != n || eigmin(A)<0) throw("RandomizedCD2 error: non-spd matrix") end

    while !OK

     i  = Int(floor(rand() * n) + 1) #rand a number between 1 and n
     Ai = view(A, i, :)

     #ei = zeros(n); ei[i] = 1.0 #unit vector in R^n
     #xk  = Ai == 0 ? x0 : x0 - dot(Ai,res)/norm(Ai,2)^2 * ei
     Aix = @kdot(n, Ai, state.x)
     if A[i,i] != 0
        state.x[i] -= (Aix - b[i])/A[i,i]
     end
     state.res .+= Ai .* ((Aix - b[i]) / A[i,i])

     OK = cheap_stop!(stp)

    end

 return stp
end
