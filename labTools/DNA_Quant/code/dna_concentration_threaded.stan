// generated with brms 2.22.0
functions {
  /* integer sequence of values
   * Args:
   *   start: starting integer
   *   end: ending integer
   * Returns:
   *   an integer sequence from start to end
   */
  array[] int sequence(int start, int end) {
    array[end - start + 1] int seq;
    for (n in 1:num_elements(seq)) {
      seq[n] = n + start - 1;
    }
    return seq;
  }
  // compute partial sums of the log-likelihood
  real partial_log_lik_lpmf(array[] int seq, int start, int end, data vector Y, data matrix Xc, vector b, real Intercept, data matrix Xc_shape, vector b_shape, real Intercept_shape, data array[] int J_1, data vector Z_1_1, vector r_1_1, data array[] int J_2, data vector Z_2_shape_1, vector r_2_shape_1) {
    real ptarget = 0;
    int N = end - start + 1;
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] shape = rep_vector(0.0, N);
    mu += Intercept + Xc[start:end] * b;
    shape += Intercept_shape + Xc_shape[start:end] * b_shape;
    for (n in 1:N) {
      // add more terms to the linear predictor
      int nn = n + start - 1;
      mu[n] += r_1_1[J_1[nn]] * Z_1_1[nn];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      int nn = n + start - 1;
      shape[n] += r_2_shape_1[J_2[nn]] * Z_2_shape_1[nn];
    }
    mu = exp(mu);
    shape = exp(shape);
    for (n in 1:N) {
      int nn = n + start - 1;
      ptarget += gamma_lpdf(Y[nn] | shape[n], shape[n] / mu[n]);
    }
    return ptarget;
  }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> Kc;  // number of population-level effects after centering
  int<lower=1> K_shape;  // number of population-level effects
  matrix[N, K_shape] X_shape;  // population-level design matrix
  int<lower=1> Kc_shape;  // number of population-level effects after centering
  int grainsize;  // grainsize for threading
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  array[N] int<lower=1> J_2;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_shape_1;
  int prior_only;  // should the likelihood be ignored?
  
  //Prior definitions #
  vector[2] beta_prior; // Students T - c(mean, sd) - df is defaulted to 3
  vector[2] betaShape_prior; // Students T - c(mean, sd) - df is defaulted to 3
  vector[2] intercept_prior; // Students T - c(mean, sd) - df is defaulted to 3
  vector[2] inteceptShape_prior; // Students T - c(mean, sd) - df is defaulted to 3
  vector[2] var_prior; // Gamma - c(a, b)
  vector[2] varShape_prior; // Gamma - c(a, b)
}
transformed data {
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  matrix[N, Kc_shape] Xc_shape;  // centered version of X_shape without an intercept
  vector[Kc_shape] means_X_shape;  // column means of X_shape before centering
  array[N] int seq = sequence(1, N);
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
  for (i in 2:K_shape) {
    means_X_shape[i - 1] = mean(X_shape[, i]);
    Xc_shape[, i - 1] = X_shape[, i] - means_X_shape[i - 1];
  }
}
parameters {
  vector[Kc] b;  // regression coefficients
  real Intercept;  // temporary intercept for centered predictors
  vector[Kc_shape] b_shape;  // regression coefficients
  real Intercept_shape;  // temporary intercept for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  array[M_1] vector[N_1] z_1;  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  array[M_2] vector[N_2] z_2;  // standardized group-level effects
}
transformed parameters {
  vector[N_1] r_1_1;  // actual group-level effects
  vector[N_2] r_2_shape_1;  // actual group-level effects
  real lprior = 0;  // prior contributions to the log posterior
  r_1_1 = (sd_1[1] * (z_1[1]));
  r_2_shape_1 = (sd_2[1] * (z_2[1]));
  lprior += student_t_lpdf(b | 3, beta_prior[1], beta_prior[2]);
  lprior += student_t_lpdf(Intercept | 3, intercept_prior[1], intercept_prior[2]);
  lprior += student_t_lpdf(b_shape | 3, betaShape_prior[1], betaShape_prior[2]);
  lprior += student_t_lpdf(Intercept_shape | 3, inteceptShape_prior[1], inteceptShape_prior[2]);
  lprior += gamma_lpdf(sd_1 | var_prior[1], var_prior[2]);
  lprior += gamma_lpdf(sd_2 | varShape_prior[1], varShape_prior[2]);
}
model {
  // likelihood including constants
  if (!prior_only) {
    target += reduce_sum(partial_log_lik_lpmf, seq, grainsize, Y, Xc, b, Intercept, Xc_shape, b_shape, Intercept_shape, J_1, Z_1_1, r_1_1, J_2, Z_2_shape_1, r_2_shape_1);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(z_2[1]);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // actual population-level intercept
  real b_shape_Intercept = Intercept_shape - dot_product(means_X_shape, b_shape);
}
