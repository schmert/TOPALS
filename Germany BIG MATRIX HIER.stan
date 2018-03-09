data {
    int A;                         // number of single-year ages
    int R;                         // number of regions
    int K;                         // number of TOPALS offsets per schedule
    
    int nlev;                     // number of levels in spatial hierarchy (should be 3)
    
    int neff[nlev];               // number of effects by level (16,115,314)

    matrix[A,K] B;                 // spline basis fns for TOPALS

    matrix[R,R-1] S;              // spatial basis fns (S %*% {iid N(0,1)}) -> CAR effects)

    matrix[A,2] lambda_star;       // standard schedules: A x 2 (male,female)

    matrix[A,R] dstar;             // expected M-F difference in log mx by (age, region)

/*--------
Stan syntax means that identically-organized arrays of 
integers and reals need different-looking declarations. 
In both cases below, varname[s][a,r] will mean
(sex s, age a, region r).  
In R, the input array in both cases is 2 x A x R
--------*/

    int D[2,A,R];                // deaths: 2-array of deaths by age,region
    matrix[A,R]  N[2];           // exposure: 2-array of person-yrs by age,region
}

parameters {
  matrix[K,R-1] eps[2];   // deep spatial params: 2-array of (R-1) x K matrices
  vector[K]     mu[2];    // global means for alphas, by component and sex: 2-array of 1xK vectors
  
  real<lower=0> sigma_car;  // sd of CAR spatial effects
  real<lower=0> sigma_sex;  // sd of sex-difference prior 
}

transformed parameters {
  matrix[K,R] alpha[2];    // TOPALS offsets: 2-array of K x R matrices
  matrix[A,R] lambda[2];   // logmx rates: 2-array of A x R matrices
  matrix[A,R] log_Dhat[2]; // log expected deaths: 2-array of A x R matrices
  
  for (s in 1:2) {
    alpha[s]  = rep_matrix(mu[s], R) + sigma_car * eps[s] * S';      // R x K 
    lambda[s] = rep_matrix(lambda_star[,s], R) + B * alpha[s];        // R x A 
    
    log_Dhat[s] = log(N[s]) + lambda[s]; 
  } 
}

model {
  matrix[A,R]   sex_diffs;
  matrix[K-1,R] alpha_diffs[2];   // 2-array of alpha(2:K)-alpha(1:(K-1)) for each region

/*----------------------
PRIORS
-----------------------*/

//--- deep parameters
  for (s in 1:2) {
    to_vector(eps[s]) ~ normal(0,1);
  }

//--- smoothing (alpha-diff) prior

  for (s in 1:2) {
    alpha_diffs[s] = alpha[s][2:K,] - alpha[s][1:(K-1),];

    to_vector(alpha_diffs[s]) ~ normal(0, 0.7071);
  }



//--- sex-difference prior

  sex_diffs = ( lambda[1] - lambda[2] ) - dstar;  

  to_vector(sex_diffs) ~ normal(0, sigma_sex);

/*----------------------
LIKELIHOOD
-----------------------*/
  for (s in 1:2) {
    to_array_1d(D[s]) ~ poisson_log( to_array_1d(log_Dhat[s]) );
  }

}