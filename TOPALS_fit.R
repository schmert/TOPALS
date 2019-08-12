#------------------------------------------------------------
# TOPALS fitting function for mortality schedules
#
# Carl Schmertmann
#   created 01 Mar 2018
#   edited  12 Aug 2019 

# + added fitted log mort, basis, etc to detailed output
# + added ability to fit age-grouped as well as single-yr data
# + changed fitting algorithm from Newton-Raphson 
#     to (penalized) IRLS
#
# Fits TOPALS parameters to single-year or age-group (D,N) 
# data by penalized IRLS with analytical derivatives
#
# A more complete explanation is in 
# https://github.com/schmert/TOPALS/blob/master/TOPALS_fitting_with_grouped_data.pdf
#------------------------------------------------------------

TOPALS_fit = function( N, D, std,
                       group_lower_age    = 0:99,
                       group_upper_age    = 1:100,
                       knot_positions     = c(0,1,10,20,40,70), 
                       penalty_precision  = 2,
                       max_iter           = 20,
                       alpha_tol          = .00005,
                       details            = FALSE) {

    require(splines)
  
    ## single years of age from 0 to (A-1)
    A   = length(std)
    age = 0:(A-1)

    ## B is an AxK matrix. Each column is a linear B-spline basis function
    B      = bs( age, knots=knot_positions, degree=1 )
    K = ncol(B) 
    
    D1 = diff( diag(K), diff=1)
    P  = penalty_precision * crossprod(D1)
    
    ## number and width of age groups
    G     = length(group_lower_age)   
    nages = group_upper_age - group_lower_age
    
    ## weighting matrix for mortality rates (assumes uniform
    ## distribution of single-year ages within groups)
    W = matrix(0, nrow=G, ncol=A, 
               dimnames=list(group_lower_age , age))

    offset = 0
    for (g in 1:G) {
      W[g, offset + 1:nages[g]] = 1/nages[g]
      offset = offset + nages[g]
    }

    ## penalized log lik function
    Q = function(alpha) {
      M = W %*% exp( std + B %*% alpha)
      likelihood = sum(D * log(M) - N * M)
      penalty    = 1/2 * t(alpha) %*% P %*% alpha
      return( likelihood - penalty )
    }
    
    #------------------------------------------------
    # iteration function: 
    # next alpha vector as a function of current alpha
    #------------------------------------------------
    next_alpha = function(alpha) {
      mu = as.vector( exp( std + B %*% alpha))
      M  = as.vector( W %*% mu)
      
      Dhat = N * M
      
      X = W %*% diag(mu) %*% B
      A = diag(N/M)
      
      y = (D-Dhat)/N + X %*% alpha
      
      updated_alpha = solve( t(X) %*% A %*% X + P, t(X) %*% A %*% y)
      return(as.vector(updated_alpha))
    }
    
    ## main iteration:     
    a = rep(0, K)
    
    niter = 0
    repeat {
      niter      = niter + 1
      last_param = a
      a          = next_alpha( a )  # update
      change     = a - last_param

      converge = all( abs(change) < alpha_tol)
      overrun  = (niter == max_iter)
      
      if (converge | overrun) { break }
      
    } # repeat
    
    if (details | !converge | overrun) {
      if (!converge) print('did not converge')
      if (overrun) print('exceeded maximum number of iterations')
      
      mu    = as.vector( exp(std + B %*% a))
      M     = as.vector( W %*% mu )
      dhat  = N * M
      
      X     = W %*% diag(mu) %*% B
      A     = diag(N/M)
      
      covar = solve( t(X) %*% A %*% X + P)
        
      return( list( alpha    = a, 
                    D        = D,
                    N        = N,
                    L        = group_lower_age,
                    U        = group_upper_age,
                    knots    = knot_positions,
                    std      = std,
                    B        = B,
                    logm     = std + B %*% a,
                    covar    = covar,
                    Qvalue   = Q(a),
                    converge = converge, 
                    maxiter  = overrun))
    } else return( a) 
    
} # TOPALS_fit

