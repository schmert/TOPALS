#------------------------------------------------------------
# TOPALS fitting function for mortality schedules
#
# Carl Schmertmann
#   created 01 Mar 2018
#   edited  16 Mar 2019 (added fitted log mort, basis, etc to detailed output)
#
# Fits TOPALS parameters to single-year (D,N) data by
# Newton-Raphson iteration with analytical derivatives
#
# A more complete explanation is in 
# https://github.com/schmert/TOPALS/blob/master/TOPALS_fitting.pdf
#------------------------------------------------------------

TOPALS_fit = function( N, D, std,
                       max_age        = 99,
                       knot_positions = c(0,1,10,20,40,70), 
                       smoothing_k    = 1,
                       max_iter       = 20,
                       alpha_tol      = .00005,
                       details        = FALSE) {

    require(splines)
  
    ## single years of age from 0 to max_age
    age = 0:max_age
    
    ## B is an Ax7 matrix. Each column is a linear B-spline basis function
    B      = splines::bs( age, knots=knot_positions, degree=1 )
    nalpha = ncol(B) 
    
    ## penalized log lik function
    Q = function(alpha) {
      lambda.hat = as.numeric( std + B %*% alpha)
      penalty    = smoothing_k * sum( diff(alpha)^2 )
      return( sum(D * lambda.hat - N * exp(lambda.hat)) - penalty)
    }
    
    ## expected deaths function
    Dhat = function(alpha) {
      lambda.hat = std + B %*% alpha
      return(  as.numeric( N * exp(lambda.hat) ))
    }      
    
    ## S matrix for penalty
    S = matrix(0,nalpha-1,nalpha) 
    diag(S[, 1:(nalpha-1)]) = -1
    diag(S[, 2:(nalpha)  ]) = +1
    SS = crossprod(S)
    
    #------------------------------------------------
    # iteration function: 
    # next alpha vector as a function of current alpha
    #------------------------------------------------
    next_alpha = function(alpha) {
      dhat = Dhat(alpha)
      M = solve ( t(B) %*% diag(dhat) %*% B + 2*smoothing_k *SS)
      v = t(B) %*% (D - dhat) - 2* (smoothing_k * (SS %*% alpha))
      return( alpha + M %*% v)
    }
    
    ## main iteration:     
    a = rep(0, nalpha)
    
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
      
      dhat = Dhat(a)
      covar = solve( t(B) %*% diag(dhat) %*% B + 2*smoothing_k *SS)
      
      return( list( alpha    = a, 
                    knots    = knot_positions,
                    std      = std,
                    B        = B,
                    offset   = B %*% a,
                    logm     = std + B %*% a,
                    covar    = covar,
                    Qvalue   = Q(a),
                    converge = converge, 
                    maxiter  = overrun))
    } else return( a) 
    
} # TOPALS_fit

