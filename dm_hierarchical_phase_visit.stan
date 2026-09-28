
data {
  int<lower=1> N_subj;
  int<lower=1> N_sess;
  int<lower=1> T;

  array[N_sess] int<lower=1, upper=N_subj> subj;
  array[N_sess] int<lower=0, upper=1> phase;

  array[N_sess, T] int<lower=0, upper=1> y;
  array[N_sess] int<lower=0,upper=1> visit2;

  array[T] real rew1;
  array[T] real rew2;
  array[T] real eff1;
  array[T] real eff2;
}

parameters {
  vector[2] mu_base;                 // baseline means: [reward, effort]
  vector[2] mu_visit;
  vector[2] mu_phase;                // phase means:    [reward, effort]

  vector<lower=0>[4] sigma_subj;     // SDs for 4 subject-level effects
  cholesky_factor_corr[4] L_subj;    // correlations between them
  matrix[4, N_subj] z_subj;          // subject latent z-scores
}

transformed parameters {
  matrix[4, N_subj] subj_eff;

  vector[N_sess] rewSens;
  vector[N_sess] effSens;

  subj_eff = diag_pre_multiply(sigma_subj, L_subj) * z_subj;

  for (n in 1:N_sess) {
    int s = subj[n];
    real ph = phase[n];

    real theta_rew =
      mu_base[1] +
      subj_eff[1, s] +
      (mu_phase[1] + subj_eff[3, s]) * ph + mu_visit[1]*visit2[n];

    real theta_eff =
      mu_base[2] +
      subj_eff[2, s] +
      (mu_phase[2] + subj_eff[4, s]) * ph + mu_visit[2]*visit2[n];

    rewSens[n] = exp(theta_rew);
    effSens[n] = exp(theta_eff);
  }
}

model {
  mu_base ~ normal(0, 1);
  mu_phase ~ normal(0, 0.5);
  mu_visit ~ normal(0,0.5);

  sigma_subj ~ exponential(1);
  L_subj ~ lkj_corr_cholesky(2.0);

  to_vector(z_subj) ~ std_normal();

  for (n in 1:N_sess) {
    for (t in 1:T) {
      real V1;
      real V2;

      V1 = rewSens[n] * rew1[t] - effSens[n] * square(eff1[t]);
      V2 = rewSens[n] * rew2[t] - effSens[n] * square(eff2[t]);

      y[n, t] ~ bernoulli_logit(V2 - V1);
    }
  }
}

generated quantities {
  real delta_rew = mu_phase[1];
  real delta_eff = mu_phase[2];

  real delta_rew_exp = exp(mu_phase[1]);
  real delta_eff_exp = exp(mu_phase[2]);

  corr_matrix[4] Omega_subj;
  vector[N_sess] rewSens_sess;
  vector[N_sess] effSens_sess;

  real ybar_rep = 0;
  real ybar_obs = 0;
  int total = 0;

  Omega_subj = multiply_lower_tri_self_transpose(L_subj);

  for (n in 1:N_sess) {
    rewSens_sess[n] = rewSens[n];
    effSens_sess[n] = effSens[n];
  }

  for (n in 1:N_sess) {
    for (t in 1:T) {
      real V1;
      real V2;
      real p;
      int y_rep;

      V1 = rewSens[n] * rew1[t] - effSens[n] * square(eff1[t]);
      V2 = rewSens[n] * rew2[t] - effSens[n] * square(eff2[t]);

      p = inv_logit(V2 - V1);
      y_rep = bernoulli_rng(p);

      ybar_rep += y_rep;
      ybar_obs += y[n, t];
      total += 1;
    }
  }

  ybar_rep /= total;
  ybar_obs /= total;
}

