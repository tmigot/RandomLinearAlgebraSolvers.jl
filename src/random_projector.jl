###############################################################################
#
# Collection of methods generating random projectors.
#
# Improvements:
#  * random_matrix_4
#  * include sparse projectors
#  * rank-1 projectors
#
###############################################################################

"""
Random projector of size k x m
with a normal distribution
"""
function random_matrix_1(k, m)
  P = zeros(k, m) #init

  #normal distribution
  P = rand(k, m)

  return 1 / sqrt(k) * P
end

"""
Random projector of size k x m
with -1 or 1 both with probability 1/2
"""
function random_matrix_2(k, m)
  P = zeros(k, m) #init

  #-1 or 1 both with probability 1/2
  P = rand([-1, 1], (k, m))

  return 1 / sqrt(k) * P
end

"""
Random projector of size k x m
with -1,0,1 respectively with probability 1/6,4/6,1/6
"""
function random_matrix_3(k, m)
  P = zeros(k, m) #init

  #-1,0,1 respectively with probability 1/6,4/6,1/6
  P = rand([-1, 0, 0, 1], (k, m))

  return 1 / sqrt(k) * P
end

"""
Random projector of size k x m
with orthogonal projection on a random k-dimensional linear subspace of R^m
"""
function random_matrix_4(k, m)
  P = zeros(k, m) #init

  #orthogonal projection on a random k-dimensional
  #linear subspace of R^m
  #P = rand(k, m)
  throw("NotImplemented")
  return 1 / sqrt(k) * P
end

"""
Check the average sparsity of random projector over N random matrices of size k*m
"""
function random_projector_sparse(projector::Function, N::Int64, k::Int64, m::Int64)
  tot = 0 #nb of non-zero elements

  for i = 1:N
    tot += SparseArrays.nnz(SparseArrays.sparse(projector(k, m))) / (k * m)
  end

  return tot / N
end

"""
Check the average rank of random projector over N random matrices of size k*m
"""
function random_projector_rank(projector::Function, N::Int64, k::Int64, m::Int64)
  tot = 0 #average rank of the projector

  for i = 1:N
    tot += rank(projector(k, m))
  end

  return tot / N
end
