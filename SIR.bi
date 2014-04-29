/**
 * SIR model.
 */
model SIR {
  const h = 0.02; // time step
  const delta_abs = 1.0e-2; // absolute error tolerance
  const delta_rel = 1.0e-5; // relative error tolerance
  
  dim n(3);
  param theta_beta[n], theta_nu[n];  // params of OU processes over fluxes
  noise w_beta, w_nu;  // noise terms
  state ln_beta, ln_nu;  // fluxes
  state s, i, r;  // susceptible, infectious, recovered
  obs y_s, y_i, y_r;  // observations
  
  inline epsilon_s = delta_abs + delta_rel*abs(s);
  inline epsilon_i = delta_abs + delta_rel*abs(i);
  inline epsilon_r = delta_abs + delta_rel*abs(r);

  sub parameter {
    // uninformative
  }

  sub proposal_parameter {
    const scale = 0.5;

    input U_beta[n,n], U_nu[n,n];  // Cholesky factors for adapted proposals
    param z_beta[n](has_output = 0), z_nu[n](has_output = 0);

    z_beta[i] ~ normal(0.0, 1.0);
    z_nu[i] ~ normal(0.0, 1.0);

    theta_beta <- theta_beta + scale*U_beta*z_beta;
    theta_nu <- theta_nu + scale*U_nu*z_nu;
  }

  sub initial {
    input s0, i0, r0;

    s <- s0;
    i <- i0;
    r <- r0;

    /* initialise rates from stationary distribution */
    ln_beta ~ normal(theta1_beta/theta2_beta, 0.5*theta3_beta**2/theta2_beta);
    ln_nu ~ normal(theta1_nu/theta2_nu, 0.5*theta3_nu**2/theta2_nu);
  }

  sub transition(delta = h) {
    w_beta ~ normal(0.0, sqrt(h));
    w_nu ~ normal(0.0, sqrt(h));
    ode(h = h, atoler = delta_abs, rtoler = delta_rel, alg = 'RK4(3)') {
      ds/dt = -exp(ln_beta)*s*i;
      di/dt = exp(ln_beta)*s*i - exp(ln_nu)*i;
      dr/dt = exp(ln_nu)*i;
      dln_beta/dt = (theta1_beta - theta2_beta*ln_beta) + theta3_beta*w_beta/h;
      dln_nu/dt = (theta1_nu - theta2_nu*ln_nu) + theta3_nu*w_nu/h;
    }
  }

  sub bridge {
    const lambda = 2.0;  // tempering

    input s_ell2, s_sf2;
    input i_ell2, i_sf2;
    input r_ell2, r_sf2;

    inline s_k = s_sf2*exp(-0.5*(t_next_obs - t_now)**2/s_ell2);
    inline s_mu = log(s)*s_k/s_sf2;
    inline s_sigma = sqrt(s_sf2 - s_k*s_k/s_sf2 + epsilon_s**2);

    y_s ~ log_normal(s_mu, lambda*s_sigma);

    inline i_k = i_sf2*exp(-0.5*(t_next_obs - t_now)**2/i_ell2);
    inline i_mu = log(i)*i_k/i_sf2;
    inline i_sigma = sqrt(i_sf2 - i_k*i_k/i_sf2 + epsilon_i**2);

    y_i ~ log_normal(i_mu, lambda*i_sigma);

    inline r_k = r_sf2*exp(-0.5*(t_next_obs - t_now)**2/r_ell2);
    inline r_mu = log(r)*r_k/r_sf2;
    inline r_sigma = sqrt(r_sf2 - r_k*r_k/r_sf2 + epsilon_r**2);

    y_r ~ log_normal(r_mu, lambda*r_sigma);
  }

  sub observation {
    y_s ~ uniform(s - epsilon_s, s + epsilon_s);
    y_i ~ uniform(i - epsilon_i, i + epsilon_i);
    y_r ~ uniform(r - epsilon_r, r + epsilon_r);
  }
}
