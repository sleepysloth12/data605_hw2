lu_decomp = function(A) {
  
  #Getting the dimensions of matrix
  
  n = dim(A)[1]
  
  #checking if matrix is square
  if (n != dim(A)[2]) {
    stop("The matrix must be square.")
  }
  
  #establishing upper and lower triangle matrix (to be filled)
  L = diag(n)
  U = matrix(0, nrow = n, ncol = n)
  
  #hard coding each LU decomp individually for 2x2 3x3 and 4x4
  
  
  switch(as.character(n),
         #2x2
         '2' = {
           L[2, 1] = A[2, 1] / A[1, 1]
           U[1,] = A[1,]
           U[2, 2] = A[2, 2] - L[2, 1] * A[1, 2]
         },
         #3x3
         '3' = {
           L[2, 1] = A[2, 1] / A[1, 1]
           L[3, 1] = A[3, 1] / A[1, 1]
           L[3, 2] = (A[3, 2] - L[3, 1] * A[1, 2]) / A[2, 2]
           
           U[1, ] = A[1, ]
           U[2, 2:3] = A[2, 2:3] - L[2, 1] * A[1, 2:3]
           U[3, 3] = A[3, 3] - L[3, 1] * A[1, 3] - L[3, 2] * A[2, 3]
           U[2, 1] = A[2, 1]
           U[3, 1:2] = A[3, 1:2] - L[3, 1] * A[1, 1:2] - L[3, 2] * A[2, 1:2]
         },
         #4x4
         '4' = {
           L[2, 1] = A[2, 1] / A[1, 1]
           L[3, 1] = A[3, 1] / A[1, 1]
           L[4, 1] = A[4, 1] / A[1, 1]
           
           L[3, 2] = (A[3, 2] - L[3, 1] * A[1, 2]) / A[2, 2]
           L[4, 2] = (A[4, 2] - L[4, 1] * A[1, 2]) / A[2, 2]
           
           L[4, 3] = (A[4, 3] - L[4, 1] * A[1, 3] - L[4, 2] * A[2, 3]) / A[3, 3]
           
           U[1,] = A[1,]
           U[2, 2:4] = A[2, 2:4] - L[2, 1] * A[1, 2:4]
           U[3, 3:4] = A[3, 3:4] - L[3, 1] * A[1, 3:4] - L[3, 2] * A[2, 3:4]
           U[4, 4] = A[4, 4] - L[4, 1] * A[1, 4] - L[4, 2] * A[2, 4] - L[4, 3] * A[3, 4]
           U[2, 1] = A[2, 1]
           U[3, 1:2] = A[3, 1:2] - L[3, 1] * A[1, 1:2] - L[3, 2] * A[2, 1:2]
           U[4, 1:3] = A[4, 1:3] - L[4, 1] * A[1, 1:3] - L[4, 2] * A[2, 1:3] - L[4, 3] * A[3, 1:3]
         },
         stop("This function only works for 2x2, 3x3, and 4x4 matrices.")
  )
  
  # calculate LXU and checking results to see if its close enough to A
  
  LU = L %*% U
  
  check_result = all.equal(A, LU, tolerance = 3)
  
  if (!isTRUE(check_result)) {
    stop(paste("Decomposition failed: ", check_result))
  }
  
  print("Matrix L:")
  print(L)
  
  print("Matrix U:")
  print(U)
  
  return(list(L = L, U = U))
  
}

#testing to see if it works

# 2x2 Test Cases
A1 = matrix(c(2, 1, 1, 2), nrow = 2)
result1 = lu_decomp(A1)

A2 = matrix(c(3, 4, 2, 1), nrow = 2)
result2 = lu_decomp(A2)

# 3x3 Test Cases
A3 = matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 10), nrow = 3)
result3 = lu_decomp(A3)

A4 = matrix(c(2, 1, 4, 3, 2, 1, 1, 2, 3), nrow = 3)
result4 = lu_decomp(A4)

# 4x4 Test Cases
A5 = matrix(c(1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6, 4, 5, 6, 7), nrow = 4)
result5 = lu_decomp(A5)

A6=matrix(c(11, 12, 13, 14, 12, 13, 14, 15, 13, 14, 15, 16, 14, 15, 16, 17), nrow = 4)
result6 = lu_decomp(A6)

# 2x2 Test Cases
A7 = matrix(c(5, 2, 3, 4), nrow = 2)
result7 = lu_decomp(A7)

A8 = matrix(c(1, 3, 2, 6), nrow = 2)
result8 = lu_decomp(A8)

# 3x3 Test Cases
A9 = matrix(c(9, 6, 3, 4, 7, 1, 2, 1, 8), nrow = 3)
result9 = lu_decomp(A9)

A10 = matrix(c(2, 4, 1, 7, 5, 3, 1, 1, 1), nrow = 3)
result10 = lu_decomp(A10)

# 4x4 Test Cases
A11 = matrix(c(1, 2, 2, 3, 4, 4, 3, 2, 3, 4, 4, 1, 1, 2, 3, 4), nrow = 4)
result11 = lu_decomp(A11)

A12 = matrix(c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4), nrow = 4)
result12 = lu_decomp(A12)

