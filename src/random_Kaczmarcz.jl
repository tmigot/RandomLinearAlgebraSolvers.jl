export RandomizedKaczmarz, RandomizedBlockKaczmarz

"""
Randomized Kaczmarz

Sect. 3.3 in Gower, R. M., & Richtárik, P. (2015).
Randomized iterative methods for linear systems.
SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.

RK takes a step in the direction of the negative stochastic gradient.
This means that it is equivalent to the SGD method.
However, the stepsize choice is very special: RK chooses the stepsize which leads to the point
which is closest to x^* in the Euclidean norm.
"""
function RandomizedKaczmarz(stp :: AbstractStopping; kwargs...)
    A,b = get_matrix(stp.pb), get_vector(stp.pb)
    m,n = size(A)
    x0  = stp.current_state.x
    xk  = x0
    stp.current_state.res = A*xk - b
    OK = start!(stp)

    while !OK

     i  = Int(floor(rand() * m)+1) #rand a number between 1 and m
     Ai = A[i,:]
     xk  = Ai == 0 ? x0 : x0 - (dot(Ai,x0)-b[i])/dot(Ai,Ai) * Ai

     update!(stp.current_state, x = xk)
     stp.current_state.res = A*xk - b
     OK = stop!(stp)
     x0  = xk

    end

 return stp
end

"""
Randomized block Kaczmarz

Sect. 3.5 in Gower, R. M., & Richtárik, P. (2015).
Randomized iterative methods for linear systems.
SIAM Journal on Matrix Analysis and Applications, 36(4), 1660-1690.
"""
function RandomizedBlockKaczmarz(stp :: AbstractStopping;
                                 r :: Int = 15,
                                 rand_r :: Bool = false,
                                 kwargs...)
    A,b = get_matrix(stp.pb), get_vector(stp.pb)
    m,n = size(A)
    x0  = stp.current_state.x
    xk  = x0
    stp.current_state.res = A*xk - b
    OK = start!(stp)

    while !OK#norm(res) > prec

     #r   = Int(floor(rand() * min(m,n))+1) #rand a number between 1 and m
     r    = rand_r ? Int(floor(rand() * m)+1) : min(r, m)
     sub = zeros(Int64, r)
     #sample!(1:m, sub) #x is a subset of [1,...,m] of size r #sample repeat entries
     sub = StatsBase.randperm(m)[1:r]

     Ai   = A[sub,:]
     resk = stp.current_state.res[sub]
     xk   = Ai == 0 ? x0 : x0 - Ai' * pinv(Matrix(Ai*Ai')) * resk

     update!(stp.current_state, x = xk)
     stp.current_state.res = A*xk - b
     OK = stop!(stp)
     x0  = xk

    end

 return stp
end
