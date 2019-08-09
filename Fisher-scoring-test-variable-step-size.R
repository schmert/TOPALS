library(HMDHFDplus)
library(tidyverse)
library(splines)

NEED.CREDENTIALS = FALSE     # switch to FALSE after credentials entered once
if (NEED.CREDENTIALS) {
  id = userInput()
  pw = userInput()
}  


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
## for single ages 0..99: Canada Females 1959

CAN = readHMDweb('CAN', 'fltper_1x1',
                   password = pw, username=id) %>%
  filter(Year==1959, Age<100) %>%
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

##==================
## sample experiment
Q = function(alpha) {
  mu1 = exp( std + B %*% alpha )
  mu = W %*% mu1
  logL    = sum( -this_N * mu + this_D * log(mu))
  penalty = - 1/2 * t(alpha) %*% P %*% alpha
  return(unlist( list( logL=logL, penalty=penalty, obj_fn = logL + penalty)))
}

## sample experiment
Qval = function(alpha) {
  mu1 = exp( std + B %*% alpha )
  mu = W %*% mu1
  return( sum( -this_N * mu + this_D * log(mu)) - 
               1/2 * t(alpha) %*% P %*% alpha )
}




factor = 1
this_N = N * factor
this_D = rpois(length(this_N), this_N * D/N)

##==================

## first exploratory iterations
alpha = rep(0,K)

D1 = diff( diag(K), diff=1) # first differencing matrix
P  = 2 * crossprod(D1)    # roughness penalty is alpha' [P] alpha

next_alpha = function(alpha) {
  mu1 = as.vector( exp( std + B %*% alpha))   # single-year rates
  mu  = as.vector( W %*% mu1)                 # age group avg rates
  
  BMW  = t(B) %*% diag(mu1) %*% t(W) %*% diag(1/mu)
  Dhat = this_N * mu 
  
  # score (gradient of penalized likelihood wrt TOPALS alphas)
  S = ( BMW %*% (this_D - Dhat) ) - P %*% alpha      
  
  # expected value of Hessian at current alpha values
  H = +( BMW %*% diag(this_N*mu) %*% t(BMW) + tcrossprod(P %*% alpha) )
  
  update = solve(H) %*% S

  ## try several step sizes
  step_size = 2^(-(0:4))
  trial_Q   = NA*step_size
  for (i in seq(step_size)) {
    trial_Q[i] = Qval( alpha + step_size[i] * update)
  }
  ibest = which.min(trial_Q)
  best_update = update * step_size[ibest]  

  new_value = alpha + best_update
  
  print(list(stepsize=step_size[ibest],
             params = data.frame(alpha, best_update, new_value),
             Qval   = data.frame(Q.old=Q(alpha), Q.new=Q(new_value))))
  
  return(new_value)
}

maxiter = 200
a = matrix(NA, K, maxiter)
a[,1] = alpha

i = 2
while (i <= maxiter) {
  a[,i] = next_alpha( a[,i-1])
  delta_a = a[,i] - a[,i-1]
  if (all(abs(delta_a) < .00005)) break
  i = i+1
}

if (i>maxiter) {
   print('--- DID NOT CONVERGE ---')
   alpha_hat = a[,maxiter]
} else {  
  alpha_hat = a[,i]
}  


## plot data

hues = c('red','darkgreen','royalblue','orangered','salmon','lawngreen')

plot(  age+.50, rate1$logmx[1:100], pch=16, ylim=c(-10,0),
       main=paste(sum(this_D),'deaths among', round(sum(this_N)),'women'),
       sub = ifelse(i>maxiter,'DID NOT CONVERGE','CONVERGED'))
lines( age+.50, std, col='grey', lwd=3)
this_color = sample(hues,1)
lines( age+.50, std + B %*% alpha_hat, col=this_color, lwd=5)

for (i in seq(L)) {
   y = log(this_D[i]/this_N[i])
   H = ifelse(i<length(L), L[i+1], 100)
   segments(L[i],y,H,y, col=this_color, lwd=4)
}

abline(v=c(L,100),lty=2, col='grey')
text(L, -9, this_D, cex=.80)
points( knot_positions, rep(-10,length(knot_positions)), pch=15, col=2,cex=1.2)

