data {
    int A;                         // number of single-year ages
    int R;                         // number of regions
    int K;                         // number of TOPALS offsets per schedule
    row_vector[A] lambda_star[2];  // standard schedule: 2-array of row vectors

    matrix[3,A] dstar;             // expected M-F difference in log mx by (zone, age) 
    int[R] zone;                   // R-vector indicating the 'zone' (E=1,W=2,Berlin=3) for each region 
}
parameters {
  matrix[R-1,K] eps[2];   // deep spatial params: 2-array of (R-1) x K matrices
  row_vector[K] mu[2];    // global means for alphas, by component and sex: 2-array of 1xK vectors

  real<lower=0> sigma_car;  // sd of CAR spatial effects
  real<lower=0> sigma_sex;  // sd of sex-difference prior 
}

transformed parameters {
  matrix[R,K] alpha[2];   // TOPALS offsets: 2-array of R x K matrices
  matrix[R,A] lambda[2];  // logmx rates: 2-array of R x A matrices

  for (s in 1:2) {
     alpha[s]  = to_matrix(mu[s], R, K, 0) + sigma_car * Ztilde * eps[s];  // R x K matrix of offsets
     lambda[s] = to_matrix(lamdbda_star[s], R, A, 0) + alpha[s] * B';      // R x A matrix of log mx
  }
}

model {
  matrix[R,A]   sex_diffs;
  matrix[R,K-1] alpha_diffs[2];   // 2-array of alpha(2:K)-alpha(1:(K-1)) for each region

/*----------------------
        PRIORS
-----------------------*/

//--- deep parameters
  to_array_1d(eps) ~ normal(0,1);

//--- smoothing (alpha-diff) prior

  for (s in 1:2) {
     alpha_diffs[s] = alpha[s][, 2:K] - alpha[s][, 1:(K-1)];
  }

  to_array_1d(alpha_diffs) ~ normal(0, 0.7071);


//--- sex-difference prior

  sex_diffs = ( lambda[1] - lambda[2] ) - dstar[ zone , ];  // @@ probably have to work on this

  to_vector(sex_diffs) ~ normal(0, sigma_sex);

/*----------------------
       LIKELIHOOD
-----------------------*/

D ~ ... 

}