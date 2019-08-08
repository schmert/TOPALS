library(HMDHFDplus)
library(tidyverse)
library(splines)

#id = userInput()
#pw = userInput()

## get Italy 1980 Female data: deaths, expos, mx for 
## both 1- and 5-year groups

rate1 = readHMDweb('ITA', 'fltper_1x1',
                  password = pw, username=id) %>%
          filter(Year==1980) %>%
          select(age=Age,mx=mx) %>%
          mutate(logmx = log(mx))

rate5 = readHMDweb('ITA', 'fltper_5x1',
                   password = pw, username=id) %>%
          filter(Year==1980) %>%
          select(age=Age,mx=mx) %>%
          mutate(logmx = log(mx))

N5 = readHMDweb('ITA', 'Exposures_5x1',
                   password = pw, username=id) %>%
  filter(Year==1980) %>%
  select(age=Age,expos=Female)

D5 = readHMDweb('ITA', 'Deaths_5x1',
                password = pw, username=id) %>%
  filter(Year==1980) %>%
  select(age=Age,deaths=Female)

ITA = full_join(N5,D5, by='age') %>%
        full_join(rate5, by='age') %>%
      filter(age < 100)

## get an arbitrary standard 
## for single ages 0..99: Canada Females 2000

CAN = readHMDweb('CAN', 'fltper_1x1',
                   password = pw, username=id) %>%
  filter(Year==2000, Age<100) %>%
  select(age=Age,mx=mx) %>%
  mutate(logmx = log(mx))

std = smooth.spline(CAN$logmx, nknots=20)$y

##################################################

age = 0:99
knot_positions = c(0,1,10,20,40,70)

## B is an Ax7 matrix. Each column is a linear B-spline basis function
B      = bs( age, knots=knot_positions, degree=1 )
K = ncol(B) 

####################################
L = c(0,1,seq(5,95,5))

W = matrix(0, nrow=length(L), ncol=length(age),
           dimnames=list(L,age))

W['0','0'] = 1
W['1',2:5] = 0.25
for (i in 3:21) W[i, 5*(i-2)+1:5] = 0.2

N = ITA$expos
D = ITA$deaths

## first exploratory iterations
alpha = rep(0,K)

eta = as.vector( exp( std + B %*% alpha))
mu  = as.vector( W %*% eta)

Dhat = N * mu

X = W %*% diag(eta) %*% B
A = diag(N/mu)

y = (D-Dhat)/N + X %*% alpha

alpha2 = solve( t(X) %*% A %*% X, t(X) %*% A %*% y)

## do it again

eta = as.vector( exp( std + B %*% alpha2))
mu  = as.vector( W %*% eta)

Dhat = N * mu

X = W %*% diag(eta) %*% B
A = diag(N/mu)

y = (D-Dhat)/N + X %*% alpha2

alpha3 = solve( t(X) %*% A %*% X, t(X) %*% A %*% y)

## do it again

eta = as.vector( exp( std + B %*% alpha3))
mu  = as.vector( W %*% eta)

Dhat = N * mu

X = W %*% diag(eta) %*% B
A = diag(N/mu)

y = (D-Dhat)/N + X %*% alpha3

alpha4 = solve( t(X) %*% A %*% X, t(X) %*% A %*% y)


## do it again

eta = as.vector( exp( std + B %*% alpha4))
mu  = as.vector( W %*% eta)

Dhat = N * mu

X = W %*% diag(eta) %*% B
A = diag(N/mu)

y = (D-Dhat)/N + X %*% alpha4

alpha5 = solve( t(X) %*% A %*% X, t(X) %*% A %*% y)

## do it again

eta = as.vector( exp( std + B %*% alpha5))
mu  = as.vector( W %*% eta)

Dhat = N * mu

X = W %*% diag(eta) %*% B
A = diag(N/mu)

y = (D-Dhat)/N + X %*% alpha5

alpha6 = solve( t(X) %*% A %*% X, t(X) %*% A %*% y)


## do it again

eta = as.vector( exp( std + B %*% alpha6))
mu  = as.vector( W %*% eta)

Dhat = N * mu

X = W %*% diag(eta) %*% B
A = diag(N/mu)

y = (D-Dhat)/N + X %*% alpha6

alpha7 = solve( t(X) %*% A %*% X, t(X) %*% A %*% y)



if (FALSE) 
{TOPALS_fit = function( N, D, std,
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

} # if FALSE