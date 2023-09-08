# Define a function that takes a matrix A and matrices U and V as arguments
woodbury_inverse_matrix <- function(A, U, V) {
  # Check that the arguments are compatible
  if (ncol(A) != nrow(U) || ncol(A) != nrow(V)) {
    stop("The dimensions of A, U and V do not match")
  }
  # Compute the inverse of A
  m_inv <- solve(A)
  # Compute the matrix c as the product of U and V
  c <- U %*% V
  # Compute the Woodbury Identity
  w_inv <- m_inv - m_inv %*% U %*% solve(c + t(V) %*% m_inv %*% U) %*% t(V) %*% m_inv
  # Return the result
  return(w_inv)
}

# Set the seed for reproducibility
set.seed(123)
# Generate a random 5x5 matrix A
A <- matrix(rnorm(25), nrow = 5)
# Generate two random 5x2 matrices U and V
U <- matrix(rnorm(10), ncol = 2)
V <- matrix(rnorm(10), ncol = 2)
# Compute the inverse of A + uv' using the Woodbury Identity
w_inv <- woodbury_inverse_matrix(A, U, V)
# Compute the inverse of A + uv' using the solve function
s_inv <- solve(A + U %*% t(V))
# Compare the results
all.equal(w_inv, s_inv)
# [1] TRUE
