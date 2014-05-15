/**
 * SIR model.
 */
model SIR {
  const h = 0.02; // time step
  const delta_abs = 1.0e-2; // absolute error tolerance
  const delta_rel = 1.0e-5; // relative error tolerance
  const epsilon = 2.0e-2; // observation error tolerance
  
  dim f(2);  // no. fluxes
  dim n(3);  // no. params in each flux (OU) process

  param theta[f,n];
  noise w[f];  // noise terms
  state ln_beta[f];  // fluxes
  state s, i, r;  // susceptible, infectious, recovered
  obs y_s, y_i, y_r;  // observations
  
  sub parameter {
    theta[f,0] ~ uniform(-100.0, 100.0);
    theta[f,1] ~ gamma(2.0, 1.0);
    theta[f,2] ~ uniform(0.0, 100.0);
  }

  sub proposal_parameter {
    const scale = 0.25;

    input C[f,n,n];  // adapted covariances for each flux
    param theta0[f,n](has_output = 0);
    param C0[f,n,n](has_output = 0);

    theta0 <- theta;
    C0 <- scale*C;

    theta[f,1] ~ truncated_normal(theta0[f,1], sqrt(C0[f,1,1]), 0.0);
    theta[f,0] ~ truncated_normal(theta0[f,0] + (C0[f,1,0]/C0[f,1,1])*(theta[f,1] - theta0[f,1]), sqrt(C0[f,0,0] - C0[f,1,0]**2/C0[f,1,1]), -100.0, 100.0);
    theta[f,2] ~ truncated_normal(theta0[f,2] + (C0[f,1,2]/C0[f,1,1])*(theta[f,1] - theta0[f,1]), sqrt(C0[f,2,2] - C0[f,1,2]**2/C0[f,1,1]), 0.0, 100.0);
  }

  sub initial {
    input s0, i0, r0;

    s <- s0;
    i <- i0;
    r <- r0;

    /* initialise rates from stationary distribution */
    ln_beta[f] ~ normal(theta[f,0]/theta[f,1], sqrt(0.5*theta[f,2]**2/theta[f,1]));
  }

  sub transition(delta = h) {
    w[f] ~ normal(0.0, sqrt(h));
    ode(h = h, atoler = delta_abs, rtoler = delta_rel, alg = 'RK4(3)') {
      ds/dt = -exp(ln_beta[0])*s*i;
      di/dt = exp(ln_beta[0])*s*i - exp(ln_beta[1])*i;
      dr/dt = exp(ln_beta[1])*i;
      dln_beta[f]/dt = theta[f,0] - theta[f,1]*ln_beta[f] + theta[f,2]*w[f]/h;
    }
  }

  sub bridge {
    const lambda = 2.0;  // tempering

    input s_ell2, s_sf2;
    input i_ell2, i_sf2;
    input r_ell2, r_sf2;

    inline s_k = s_sf2*exp(-0.5*(t_next_obs - t_now)**2/s_ell2);
    inline s_mu = log(s)*s_k/s_sf2;
    inline s_sigma = sqrt(s_sf2 - s_k*s_k/s_sf2 + epsilon**2);

    y_s ~ log_normal(s_mu, lambda*s_sigma);

    inline i_k = i_sf2*exp(-0.5*(t_next_obs - t_now)**2/i_ell2);
    inline i_mu = log(i)*i_k/i_sf2;
    inline i_sigma = sqrt(i_sf2 - i_k*i_k/i_sf2 + epsilon**2);

    y_i ~ log_normal(i_mu, lambda*i_sigma);

    inline r_k = r_sf2*exp(-0.5*(t_next_obs - t_now)**2/r_ell2);
    inline r_mu = log(r)*r_k/r_sf2;
    inline r_sigma = sqrt(r_sf2 - r_k*r_k/r_sf2 + epsilon**2);

    y_r ~ log_normal(r_mu, lambda*r_sigma);
  }

  sub observation {
    y_s ~ uniform(s - epsilon, s + epsilon);
    y_i ~ uniform(i - epsilon, i + epsilon);
    y_r ~ uniform(r - epsilon, r + epsilon);
  }
}
