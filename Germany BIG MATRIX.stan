data {
    int A;                         // number of single-year ages
    int R;                         // number of regions
    int K;                         // number of TOPALS offsets per schedule
    row_vector[A] lambda_star[2];  // standard schedule: 2-array of row vectors

    matrix[A,3] dstar;             // expected M-F difference in log mx by (zone, age) 
    int[R] zone;                   // R-vector indicating the 'zone' (E=1,W=2,Berlin=3) for each region 
    
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
     alpha[s]  = to_matrix(mu[s], K, R) + sigma_car * eps[s] * Ztilde';  // R x K 
     lambda[s] = to_matrix(lamdbda_star[s], A, R) + B * alpha[s];        // R x A 
     
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
  to_array_1d(eps) ~ normal(0,1);

//--- smoothing (alpha-diff) prior

  for (s in 1:2) {
     alpha_diffs[s] = alpha[s][2:K,] - alpha[s][1:(K-1),];
  }

  to_array_1d(alpha_diffs) ~ normal(0, 0.7071);


//--- sex-difference prior

  sex_diffs = ( lambda[1] - lambda[2] ) - dstar[ , zone ];  // @@ probably have to work on this

  to_vector(sex_diffs) ~ normal(0, sigma_sex);

/*----------------------
       LIKELIHOOD
-----------------------*/

  to_array_1d(D) ~ poisson_log( log_Dhat);

}