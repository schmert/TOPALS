library(HMDHFDplus)
library(tidyverse)
library(splines)

#id = userInput()
#pw = userInput()

## if necessary, get Italy 1980 Female data: 
## deaths, expos, mx for 
## both 1- and 5-year groups

NEED.HMD.DOWNLOAD = !exists('ITA.Rdata')

if (NEED.HMD.DOWNLOAD) {
  
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
    
    save(rate1, rate5, N5, D5, ITA, std, file='ITA.Rdata')
} else {
    load('ITA.Rdata')
}
##################################################

age = 0:99
knot_positions = c(0,1,10,20,40,70)

## B is an Ax7 matrix. Each column is a linear B-spline basis function
B      = bs( age, knots=knot_positions, degree=1 )
K = ncol(B) 

D1 = diff( diag(K), diff=1)
P  = 2 * crossprod(D1)
####################################
L = c(0,1,seq(5,95,5))

W = matrix(0, nrow=length(L), ncol=length(age),
           dimnames=list(L,age))

W['0','0'] = 1
W['1',2:5] = 0.25
for (i in 3:21) W[i, 5*(i-2)+1:5] = 0.2

bigN = ITA$expos
bigD = ITA$deaths

#######################################################
## experimental sample
#######################################################
target_pop = 100
factor     = target_pop / sum(bigN)

N = bigN * factor
D = rpois(length(N), N * bigD/bigN)

Q = function(alpha) {
  mu1 = exp( std + B %*% alpha )
  mu = W %*% mu1
  logL    = sum( -N * mu + D * log(mu))
  penalty = - 1/2 * t(alpha) %*% P %*% alpha
  return(unlist( list( logL=logL, penalty=penalty, obj_fn = logL + penalty)))
}

## first exploratory iterations
alpha = rep(0,K)

next_alpha = function(alpha) {
  eta = as.vector( exp( std + B %*% alpha))
  mu  = as.vector( W %*% eta)
  
  Dhat = N * mu
  
  X = W %*% diag(eta) %*% B
  A = diag(N/mu)
  
  y = (D-Dhat)/N + X %*% alpha
  
  updated_alpha = solve( t(X) %*% A %*% X + P, t(X) %*% A %*% y)
  return(updated_alpha)
}

maxiter = 25
a = matrix(NA, K, maxiter)
a[,1] = alpha

i = 2
while (i <= maxiter) {
  a[,i] = next_alpha( a[,i-1])
  delta_a = a[,i] - a[,i-1]
  if (all(abs(delta_a) < .00005)) break
  i = i+1
}

alpha_hat = a[,i]

## plot data
op = par(no.readonly = TRUE)
par(mfrow=c(1,2))

hues = c('red','darkgreen','gold','violet','coral','royalblue','orangered','salmon','lawngreen')
this_color = sample(hues,1)

plot( age+.50, rate1$logmx[1:100], pch=16, ylim=c(-10,0),
      main=paste(sum(D),'deaths among', round(sum(N)),'women',
                 '\n',
           ifelse(i<maxiter, 'CONVERGED', 'DID NOT CONVERGE')))

lines( age+.50, std, col='grey', lwd=3)
lines( age+.50, std + B %*% alpha_hat, col=this_color, lwd=3)

for (i in seq(L)) {
   y = log(D[i]/N[i])
   H = ifelse(i<length(L), L[i+1], 100)
   segments(L[i],y,H,y, col=this_color, lwd=3)
}
abline(v=c(L,100),lty=2, col='grey')

for (i in seq(L)) {
  y = log(D[i]/N[i])
  H = ifelse(i<length(L), L[i+1], 100)
  segments(L[i],y,H,y, col=this_color, lwd=4)
}

abline(v=c(L,100),lty=2, col='grey')
text(L, c(-8.5,rep(-9,length(L)-1)), D, cex=.50)
points( knot_positions, rep(-10,length(knot_positions)), pch=15, col=2,cex=1.2)

apply(a,2,Q)

matplot(t(a), type='o', main='Offsets', xlab="Iteration")

par(op)

