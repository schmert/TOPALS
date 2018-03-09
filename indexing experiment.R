### experiment with indexing order when passing data from R -> Stan 

rm(list=ls())
library(rstan)

integer_D = array(NA, c(2,8,3))  # data
real_D    = array(NA, c(2,8,3))  # data

for (i in 1:2) {
  for (j in 1:8) {
    for (k in 1:3) {
      real_D[i,j,k] = 100*i + 10*j + k
      integer_D[i,j,k] = 100*i + 10*j + k
    }
  }
}

stanDataList = list(
  integer_D = integer_D,
  real_D    = real_D
)

stanModelText = '
data {
  int integer_D[2,8,3];
  matrix[8,3] real_D[2];
}
transformed data {
  for (i in 1:8) {
    for (j in 1:3) {
      for (k in 1:2) {
         print(\"   integer_D[\", k, \"][\", i,\" ,\", j, \"]\",  integer_D[k][i,j]);
         print(\"real_D[\", k, \"][\", i,\" ,\", j, \"]\",  real_D[k][i,j]);
      }
    }
  }
}
parameters {
  real z;
}
model {
  z ~ normal(0,1);
}
'

myModel = stan_model( model_name = 'xx', model_code = stanModelText )

fit = sampling(myModel, data=stanDataList, chains=1, iter=100)

