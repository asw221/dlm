
// Fit the model "fit2" from dlmBE's README using MC Stan
//
data {
  int<lower=0> N;        // total number of observations
  int<lower=0> P;        // total number of regression coefficients
  int<lower=0> P_cr;     // number of "random effects" cr term
  int<lower=0> P_crfem;  // number of "random effects" cr:female term
  // vector[N] y;           // outcome variable
  real y[N];
  matrix[N,P] X;         // fixed effects design matrix
}

parameters {
  vector[P - P_cr - P_crfem] beta;
  vector[P_cr] theta_cr;
  vector[P_crfem] theta_crfem;
  real<lower=0,upper=100> sigma;        // residual SD
  real<lower=0,upper=100> sigma_cr;     // RE cr SD
  real<lower=0,upper=100> sigma_crfem;  // RE cr:female SD
}

transformed parameters {
  vector[P] theta = append_row(append_row(beta, theta_cr), theta_crfem);
  vector[N] mu = X * theta;
}

model {
  beta ~ normal(0, 100);
  theta_cr ~ normal(0, sigma_cr);
  theta_crfem ~ normal(0, sigma_crfem);
  y ~ normal(mu, sigma);
}
