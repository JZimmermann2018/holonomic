library(MASS)


submultisets <- function(n, k) {

   return(choose(n+k-1, k))
}

allsubmultisets <- function(n, k) {

   if (k == 0) return(1)

   allsub = 0
   for (i in 1:k) allsub = allsub + submultisets(n, i)

   return(allsub)
}

# DEPRECATED: recommended only as a reference method, for best
#             results use NeumaierSum.
#
# Kahan summation algorithm for adding up numbers having vastly
# differing magnitudes
#
# x = vector of numbers
KahanSum <- function(x) {

   n = length(x)

   sum = 0
   c = 0
   for (i in 1:n) {
      y = x[i] - c
      t = sum + y
      c = (t - sum) - y
      sum = t
   }

   return(sum)
}

# Neumaier summation algorithm for adding up numbers having vastly
# differing magnitudes
#
# x = vector of numbers
NeumaierSum <- function(x) {

   n = length(x)

   sum = 0
   c = 0
   for (i in 1:n) {
      t = sum + x[i]
      if (abs(sum) >= abs(x[i])) {
         c = c + ((sum - t) + x[i])
      }
      else {
         c = c + ((x[i] - t) + sum)
      }
      sum = t
   }

   return(sum + c)
}


# vector field to gradient
vf2grad <- function(vf, location) {

   return(vf(location))
}

# conmpute trajectory for n steps wrt. vector field 'vr',
# initial value 'initial' and step size 'h'
trajectory <- function(vf, initial, h, n) {

   if (n==0) return(initial)
   
   trajectory = initial
   current = initial
   for (i in 1:n) {
      next_loc = current + vf(current)*h
      trajectory = cbind(trajectory, next_loc)
      current = next_loc
   }

   return(trajectory)
}


vf1 <- function(x) {

   res1 = 2*x[1] + 3*x[2]
   res2 = -x[1] + x[2]

   return(c(res1, res2))
}

vf2 <- function(x) {

   return(x^2)
}

lorenz_vf <- function(x, beta = 8/3, rho = 28, sigma = 10) {

   Dx = numeric(3)

   Dx[1] = sigma*(x[2] - x[1])
   Dx[2] = x[1]*(rho - x[3]) - x[2]
   Dx[3] = x[1]*x[2] - beta*x[3]

   return(Dx)
}

pendulum_vf <- function(x, k = .1) {

   Dx = numeric(2)

   Dx[1] = x[2]
   Dx[2] = -k*sin(x[1])

   return(Dx)
}

particle_vf <- function(x, g = 9.81) {
   
   Dx = numeric(4)

   Dx[1] = x[2]
   Dx[2] = 0
   Dx[3] = x[4]
   Dx[4] = -g

   return(Dx)
}

van_der_pol_vf <- function(x, mu=5) {

   Dx = numeric(2)

   Dx[1] = x[2]
   Dx[2] = mu*(1-x[1]^2)*x[2] - x[1]

   return(Dx)
}

polyreg <- function(X, y) {

   return(solve(t(X)%*%X)%*%t(X)%*%y)
}

Delta <- function(X) {
   m = dim(X)[1]
   n = dim(X)[2]

   if (n < 2) {
      print("ERROR: Difference operator not applicable: data series consists of less than 2 entries.")
      return()
   }

   for (i in 1:(n-1)) {
      X[,i] = X[,i+1] - X[,i]
   }

   return(X[,-n])
}

# Right shift operator
#
# Add lag vectors up to lag k to data matrix X
# wrt. to a set of trajectories
add_lag_dimensions <- function(X, k, traj_set = 1:dim(X)[1]) {
   n = dim(X)[2]
   
   if (n < k+1) {
      print(paste("ERROR: Right shift operator not applicable: data series consists of less than", k+1, "entries."))
      return()
   }

   if (k < 1) return(X)

   # case length(traj_set) == 1 needs to be handled separately because selecting
   # from a 1xn matrix results in a type cast to a vector
   if (length(traj_set) == 1) {
      Xlag = lag_dimensions(t(X[traj_set,]), k)
   }
   else {
      Xlag = lag_dimensions(X[traj_set,], k)
   }

   return(rbind(X[, -(1:k)], Xlag))
}

# Right shift operator
#
# Compute lag vectors up to lag k of data matrix X
lag_dimensions <- function(X, k) {
   n = dim(X)[2]

   if (n < k+1) {
      print(paste("ERROR: Right shift operator not applicable: data series consists of less than", k+1, "entries."))
      return()
   }

   if (k < 1) return(c())

   Xlag = c()

   if (k > 1) {
      for (i in 1:(k-1)) {
         Xlag = rbind(Xlag, X[, -c(1:(k-i), (n-i+1):n)])
      }
   }

   Xlag = rbind(Xlag, X[, -((n-k+1):n)])

   return(Xlag)
}


# DMD using pseudo inverse
DMDpi <- function(X, observables = 1:dim(X)[1]) {
   m = dim(X)[1]
   n = dim(X)[2]

   if (n < 2) {
      print("DMDpi: less than 2 data points.")
      return(diag(0))
   }

   # needs to be handled separately because selecting
   # from a 1xn matrix results in a type cast to a vector
   if (m == 1) return(X[,-1]%*%ginv(t(X[,-n])))

   return(X[observables,-1]%*%ginv(X[,-n]))
}

# Reconstruction error wrt. a set of trajectories
recoError <- function(X, A, traj_set = 1:dim(X)[1]) {

#   print(traj_set)
   m = dim(X)[1]
   n = dim(X)[2]
   mp = length(traj_set)

   E = X[,-1] - A%*%X[,-n]
   subE = as.matrix(E[traj_set,])

   return(normFrobNorm(subE))
}


# Frobenius Norm
FrobNorm <- function(A) {

   return(sqrt(sum(abs(A)^2)))
}

# normalized Frobenius Norm
normFrobNorm <- function(A) {
   m = dim(A)[1]
#   print(m)
   n = dim(A)[2]
#   print(n)

   return(FrobNorm(A)/sqrt(m*n))
}

# X is a data matrix (time is column index)
#
# 'exponents' defines a monomial over the data dimensions,
# given as a list of pairs (variable index, exponent)
#
# returns X + monomial dimension
add_monomial_dimension <- function(X, exponents) {

   m = dim(X)[1]
   n = dim(X)[2]
   d = length(exponents)

   extra_dim = X[exponents[[1]][1],]^exponents[[1]][2]

   if (d > 1) {
      for (i in 2:d) {
         extra_dim = extra_dim*X[exponents[[i]][1],]^exponents[[i]][2]
      }
   }

   return(rbind(X, extra_dim))
}


# X is a data matrix (time is column index)
#
# 'exponents' defines a monomial over the data dimensions,
# given as a list of pairs (variable index, exponent)
#
# returns the monomial dimension
monomial_dimension <- function(X, exponents) {

   m = dim(X)[1]
   n = dim(X)[2]
   d = length(exponents)

   extra_dim = X[exponents[[1]][1],]^exponents[[1]][2]

   if (d > 1) {
      for (i in 2:d) {
         extra_dim = extra_dim*X[exponents[[i]][1],]^exponents[[i]][2]
      }
   }

   return(extra_dim)
}


# Add all monomial dimensions up to a certain degree
# wrt. a set of trajectories
add_all_monomial_dimensions_up_to_degree_d <- function(X, d, traj_set = 1:dim(X)[1]) {

   if (d < 2) return(X)

   Xmonomial = all_monomial_dimensions_up_to_degree_d(X[traj_set,], d)

   return(rbind(X, Xmonomial))
}

# Compute all monomial dimensions up to a certain degree
# wrt. a set of trajectories
all_monomial_dimensions_up_to_degree_d <- function(X, d) {

   if (d < 2) return(c())

   Xmonomial = c()
   for (i in 2:d) {
      Xmonomial = rbind(Xmonomial, all_monomial_dimensions_of_degree_d(X, i))
   }

   return(Xmonomial)
}

all_monomial_dimensions_of_degree_d <- function(X, d) {

   nvar = dim(X)[1]

   mset = initial_mset(nvar, d)

   M = c()
   
   repeat {
      exponents = vector2pair_list_mset(mset)

      M = rbind(M, monomial_dimension(X, exponents))

      if (!has_next_mset(mset)) break

      mset = next_mset(mset)
   }

   return(M)
}

# generate trajectory with DMD matrix
genDMD <- function(A, initial, t) {
   m = dim(A)[1]
   n = dim(A)[2]
   
   trajectory = initial
   current = initial

   new = numeric(m)
   for (i in 1:t) {
      new = A%*%current

      trajectory = cbind(trajectory, new)
      current = new
   }

   return(trajectory)
}

# generate van der Pol trajectory with DMD matrix based on delay and monom dimension
genDMDvdp <- function(A, initial, t) {
   m = dim(A)[1]
   n = dim(A)[2]
   
   trajectory = initial
   current = initial

   new = numeric(m)
   for (i in 1:t) {
      new[1] = A[1,]%*%current
      new[2] = new[1] - current[1]
      new[3] = new[1]^2*new[2]

      trajectory = cbind(trajectory, new)
      current = new
   }

   return(trajectory)
}

# generate Lorenz trajectory with DMD matrix
genDMDlorenz <- function(A, initial, t) {
   m = dim(A)[1]
   n = dim(A)[2]
   
   trajectory = initial
   current = initial

   new = numeric(m)
   for (i in 1:t) {
      new[1] = A[1,]%*%current
      new[2] = A[2,]%*%current
      new[3] = A[3,]%*%current
      new[4] = new[1]*new[2]
      new[5] = new[1]*new[3]
      new[6] = new[2]*new[3]
      new[7] = new[1]^2
      new[8] = new[2]^2
      new[9] = new[3]^2

      trajectory = cbind(trajectory, new)
      current = new
   }

   return(trajectory)
}

# generate x-Lorenz trajectory with DMD Matrix based on decoupled recurrence equation
genDMDlorenz_decoupled <- function(A, initial, t) {
   m = dim(A)[1]
   n = dim(A)[2]
   
   trajectory = initial
   current = initial

   new = numeric(m)
   for (i in 1:t) {
      new[1] = A[1,]%*%current
      new[2] = current[1]
      new[3] = current[2]
      new[4] = new[3]^2*new[2]
      new[5] = new[3]*new[2]^2
      new[6] = new[2]*new[1]/new[3]
      new[7] = new[2]^2/new[3]

      trajectory = cbind(trajectory, new)
      current = new
   }

   return(trajectory)
}


vector2pair_list_mset <- function(s) {

   n = length(s)

   exponents = c()
   for (i in 1:n) {

      if (s[i] != 0) {
         exponents = c(exponents, list(c(i, s[i])))
      }
   }

   return(exponents)
}

initial_mset <- function(n, k) {

   s = numeric(n)
   s[n] = k

   return(s)
}

has_next_mset <- function(s) {

   if (all(s[-1] == 0)) {
      return(FALSE)
   }
   else {
      return(TRUE)
   }
}

# s = multiset of cardinality 'k' over base set of cardinality 'n',
#     represented by a vector of length 'n' where entries are the
#     multiplicities of base set elements
next_mset <- function(s) {

   n = length(s)

   if (n == 0) {
      print("next_mset: length of mset is zero.")
      return(-1)
   }

   if (!(all(s[] >= 0))) {
      print("next_mset: mset contains negative entries.")
      return(-1)
   }

   if (sum(s) == 0) {
      print("next_mset: mset is zero vector.")
      return(-1)
   }
   
   # position of last non-zero entry
   for (i in 1:n) {
      if (s[n+1-i] != 0) {
         zpos = n+1-i
         break
      }
   }

   if (zpos == 1) return(s)
   
   s[zpos-1] = s[zpos-1] + 1
   s[n] = s[zpos] - 1
   if (zpos != n) s[zpos] = 0

   return(s)
}


# generate the design matrix and response vector for the holonomic
# regression problem
#
# a = (a_0, .., a_k) is the vector of known first coefficients of
#     a Taylor series
# r = order of recursion
# d = degree of coefficient polynomials
# fixed2one = c(i,j), -1 =< i =< r, 0 =< j =< d, where b_i,j is the
#             regression weight set to 1 (otherwise the null vector
#             would be a solution). c(-1,j) are the coefficients of
#             the polynomial defining the inhomogeneous part of
#             the holonomic equation.
holonomic_model_inhomogeneous <- function(a, r, d, fixed2one = c(r, d)) {

   k = length(a)-1

   if (r > k) {
      print("holonomic_model: too few coefficients for recursion order.")
      return(-1)
   }

   n_vector = numeric(d+1)
   # X = design matrix
   X = matrix(nrow = k-r+1, ncol = (r+2)*(d+1))
   
   for (n in 0:(k-r)) {
      n_vector = n^(0:d)
      X[n+1,1:(d+1)] = -n_vector
      for (i in 2:(r+2)) {
         X[n+1, ((i-1)*(d+1)+1):(i*(d+1))] = a[n+i-1]*n_vector
      }
   }

   i = fixed2one[1]
   j = fixed2one[2]
   response = (i+1)*(d+1) + (j+1)
   
   y = -X[, response]

   X = X[, -response]

   return(list(X, y, fixed2one))
}

# generate the design matrix and response vector for the holonomic
# regression problem
#
# a = (a_0, .., a_k) is the vector of known first coefficients of
#     a Taylor series
# r = order of recursion
# d = degree of coefficient polynomials
# active = list of active regression weights, others are set to zero.
# fixed2one = number of regression weight in active list which is set
#             to 1 (otherwise the null vector would be a solution).
holonomic_model_active_subset <- function(a, r, d, active, fixed2one = length(active)) {

   k = length(a)-1
   # number of active wights
   naw = length(active)

   if (r > k) {
      print("holonomic_model: too few coefficients for recursion order.")
      return(-1)
   }

   # X = design matrix
   X = matrix(nrow = k-r+1, ncol = naw)
   
   for (n in 0:(k-r)) {
      for (l in 1:naw) {
         i = active[[l]][1]
         j = active[[l]][2]
         
         if (i == -1) {
            X[n+1, l] = -n^j
         }
         else {
            X[n+1, l] = a[n+i+1]*n^j
         }
      }
   }

   response = fixed2one
   
   y = -X[, response]

   X = X[, -response]

   return(list(X, y, active[[fixed2one]]))
}

# generate the design matrix and response vector for the holonomic
# regression problem
#
# a = (a_0, .., a_k) is the vector of known first coefficients of
#     a Taylor series
# r = order of recursion
# d = degree of coefficient polynomials
# fixed2one = c(i,j), 0 =< i =< r, 0 =< j =< d, where b_i,j is the
#             regression weight set to 1 (otherwise the null vector
#             would be a solution).
holonomic_model <- function(a, r, d, fixed2one = c(r, d)) {

   k = length(a)-1

   if (r > k) {
      print("holonomic_model: too few coefficients for recursion order.")
      return(-1)
   }

   n_vector = numeric(d+1)
   # X = design matrix
   X = matrix(nrow = k-r+1, ncol = (r+1)*(d+1))
   
   for (n in 0:(k-r)) {
      n_vector = n^(0:d)
      for (i in 1:(r+1)) {
         X[n+1, ((i-1)*(d+1)+1):(i*(d+1))] = a[n+i]*n_vector
      }
   }

   i = fixed2one[1]
   j = fixed2one[2]
   response = i*(d+1) + (j+1)
   
   y = -X[, response]

   X = X[, -response]

   return(list(X, y, fixed2one))
}

normRegressionError <- function(X, y, beta) {

   # sum of squares
   ss = sum((y - X%*%beta)^2)

   n = length(y)

   return(sqrt(ss/n))
}


evaluate_polynomial <- function(poly, x) {

   d = length(poly) - 1

   if (d == -1) return(0)

   val = poly[1]
   if (d > 0) {
      for (i in 1:d) {
         val = val + poly[i+1]*x^i
      }
   }

   return(val)
}

# computes a holonomic sequence
#
# plist = list of coefficient polynomials
# initial = vector of initial values
# n = length of sequence to be generated
holonomic_sequence <- function(plist, initial, n) {

   r = length(plist)-1

   if (length(initial) != r) {
      print("holonomic_sequence: number of initial value does not match order of recursion.")
      return(-1)
   }

   if (n == 0) return(c())

   if (n <= r) return(initial[1:n])
   
   a = numeric(n)
   a[1:r] = initial
   for (i in 1:(n-r)) {
      
      sum = 0
      for (j in 0:(r-1)) {
         sum = sum + evaluate_polynomial(plist[[j+1]], i-1)*a[i+j]
      }
      
      a[r+i] = -sum/evaluate_polynomial(plist[[r+1]], i-1)
   }

   return(a)
}


generating_function <- function(a, x) {

   n = length(a)-1

   return(sum(a*x^(0:n)))
}


# computes integrate(x^n * sin(omega*x), a, b)
inner_product_power_sine <- function(a, b, n, omega = 1) {

   c = numeric(n+1)

   c[n+1] = factorial(n)/omega^(n+1)
   if (n > 0) {
      for (i in 1:n) {
         c[n+1-i] = (a*omega)/i*c[n+2-i]
      }
   }
   
   Fa = 0
   for (i in 0:n) {
      Fa = Fa - c[i+1]*cos(omega*a + i*pi/2)
   }

   c[n+1] = factorial(n)/omega^(n+1)
   if (n > 0) {
      for (i in 1:n) {
         c[n+1-i] = (b*omega)/i*c[n+2-i]
      }
   }
   
   Fb = 0
   for (i in 0:n) {
      Fb = Fb - c[i+1]*cos(omega*b + i*pi/2)
   }

   return(Fb - Fa)
}

# computes integrate(x^n * sin(omega*x), a, b)
# using Neumaier summation
inner_product_power_sine_NM <- function(a, b, n, omega = 1) {

   c = numeric(n+1)
   
   c[n+1] = factorial(n)/omega^(n+1)
   if (n > 0) {
      for (i in 1:n) {
         c[n+1-i] = (a*omega)/i*c[n+2-i]
      }
   }
   
   Fas = numeric(n+1)
   for (i in 0:n) {
      Fas[i+1] = -c[i+1]*cos((omega*a + i*pi/2)%%(2*pi))
   }
   #   Fa = -NeumaierSum(Fas)
   #   cat("Fas =", Fas, "\n")

   c[n+1] = factorial(n)/omega^(n+1)
   if (n > 0) {
      for (i in 1:n) {
         c[n+1-i] = (b*omega)/i*c[n+2-i]
      }
   }
   
   Fbs = numeric(n+1)
   for (i in 0:n) {
      Fbs[i+1] = -c[i+1]*cos((omega*b + i*pi/2)%%(2*pi))
   }
   #   Fb = -NeumaierSum(Fbs)
   #   cat("Fbs =", Fbs, "\n")

   return(NeumaierSum(c(Fbs, -Fas)))
}


# TEST FUNCTION
#
# computes integrate(x^n * sin(omega*x), a, b)
# using Neumaier summation
#
# even entries are set to 0 (to compute exact values for
# interval (-pi, pi) and omega = 1)
inner_product_power_sine_NM_test <- function(a, b, n, omega = 1) {

   if (n%%2 == 0) return(0)

   c = numeric(n+1)
   
   c[n+1] = factorial(n)/omega^(n+1)
   if (n > 0) {
      for (i in 1:n) {
         c[n+1-i] = (a*omega)/i*c[n+2-i]
      }
   }
   
   Fas = numeric(n+1)
   for (i in 0:n) {
      Fas[i+1] = -c[i+1]*cos((omega*a + i*pi/2)%%(2*pi))
   }
   #   Fa = -NeumaierSum(Fas)
   #   cat("Fas =", Fas, "\n")

   c[n+1] = factorial(n)/omega^(n+1)
   if (n > 0) {
      for (i in 1:n) {
         c[n+1-i] = (b*omega)/i*c[n+2-i]
      }
   }
   
   Fbs = numeric(n+1)
   for (i in 0:n) {
      Fbs[i+1] = -c[i+1]*cos((omega*b + i*pi/2)%%(2*pi))
   }
   #   Fb = -NeumaierSum(Fbs)
   #   cat("Fbs =", Fbs, "\n")

   return(NeumaierSum(c(Fbs, -Fas)))
}


design_matrix_power_base <- function(a, b, n) {

   X = matrix(0, nrow = n+1, ncol = n+1)

   for (i in 0:n) {
      for (j in 0:n) {
         X[i+1, j+1] = (b^(i+j+1) - a^(i+j+1))/(i+j+1)
      }
   }

   return(X)
}

response_vector_power_sine <- function(a, b, n, omega = 1) {

   y = numeric(n+1)

   for (i in 0:n) {
      y[i+1] = inner_product_power_sine_NM(a, b, i, omega)
   }

   return(y)
}

response_vector_power_sine_test <- function(a, b, n, omega = 1) {

   y = numeric(n+1)

   for (i in 0:n) {
      y[i+1] = inner_product_power_sine_NM_test(a, b, i, omega)
   }

   return(y)
}
