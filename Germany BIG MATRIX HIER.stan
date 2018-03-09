data {
    int A;                         // number of single-year ages
    int R;                         // number of regions
    int K;                         // number of TOPALS offsets per schedule

    matrix[A,K] B;                 // spline basis fns for TOPALS
    
    int neff[3];                   // number of effects at each hier. level (16,116,314)
    
/*-------------------------
spatial effects matrices 
based on regional hierarchy
--------------------------*/
    matrix[R,  neff[1]-1 ] S1;
    matrix[R,  neff[2]-1 ] S2;
    matrix[R,  neff[3]-1 ] S3;    

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
  matrix[K,neff[1]-1] eps1[2];   // deep spatial params: 2-array of K x (neff1-1)
  matrix[K,neff[2]-1] eps2[2];   // deep spatial params: 2-array of K x (neff2-1)
  matrix[K,neff[3]-1] eps3[2];   // deep spatial params: 2-array of K x (neff3-1)
    
  vector[K]     mu[2];    // global means for alphas, by comp and sex: 2-array of 1xK vectors
  
  real<lower=0> sigma_hier[3];  // sd of CAR spatial effects
  real<lower=0> sigma_sex;      // sd of sex-difference prior 
}

transformed parameters {
  matrix[K,R] alpha[2];    // TOPALS offsets: 2-array of K x R matrices
  matrix[A,R] lambda[2];   // logmx rates: 2-array of A x R matrices
  matrix[A,R] log_Dhat[2]; // log expected deaths: 2-array of A x R matrices
  
  for (s in 1:2) {
    alpha[s]  = rep_matrix(mu[s], R) 
                       + sigma_hier[1] * eps1[s] * S1' 
                       + sigma_hier[2] * eps2[s] * S2' 
                       + sigma_hier[3] * eps3[s] * S3' ;
    
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