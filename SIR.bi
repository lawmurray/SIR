/**
 * SIR model.
 */
model SIR {
  const h = 0.1; // time step
  const delta_abs = 1.0e-2; // absolute error tolerance
  const delta_rel = 1.0e-5; // relative error tolerance
  
  input s0, i0, r0;      // initial condition
  param beta, nu;        // transmission rate, recovery rate
  param sigma1, sigma2;  // diffusion rates
  state s, i, r;         // susceptible, infectious, recovered
  noise w1, w2;          // noise terms
  obs y_s, y_i, y_r;
  
  /* bridge weighting function parameters */
  input s_ell2, s_sf2;
  input i_ell2, i_sf2;
  input r_ell2, r_sf2;
  const lambda = 1.0;

  inline epsilon_s = delta_abs + delta_rel*abs(s);
  inline epsilon_i = delta_abs + delta_rel*abs(i);
  inline epsilon_r = delta_abs + delta_rel*abs(r);

  sub parameter {
    beta ~ uniform(0.0, 1.0);
    nu ~ uniform(0.0, 1.0);
    sigma1 ~ uniform(0.0, 1.0);
    sigma2 ~ uniform(0.0, 1.0);
  }

  sub proposal_parameter {
    beta ~ truncated_gaussian(beta, 2.0e-4, 0.0, 1.0);
    nu ~ truncated_gaussian(nu, 4.0e-2, 0.0, 1.0);
    sigma1 ~ truncated_gaussian(sigma1, 1.0e-3, 0.0, 1.0);
    sigma2 ~ truncated_gaussian(sigma2, 3.0e-2, 0.0, 1.0);
  }

  sub initial {
    s <- s0;
    i <- i0;
    r <- r0;
  }

  sub transition(delta = h) {
    w1 ~ gaussian(0.0, sqrt(h));
    w2 ~ gaussian(0.0, sqrt(h));
    ode(h = h, atoler = delta_abs, rtoler = delta_rel, alg = 'RK4(3)') {
      ds/dt = -beta*i*s + sigma1*s*w1/h;
      di/dt = beta*i*s - nu*i - sigma1*s*w1/h + sigma2*i*w2/h;
      dr/dt = nu*i - sigma2*i*w2/h;
    }
  }

  sub bridge {
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
    //y_s ~ gaussian(s, 0.5*epsilon_s);
    //y_i ~ gaussian(i, 0.5*epsilon_i);
    //y_r ~ gaussian(r, 0.5*epsilon_r);

    y_s ~ uniform(s - epsilon_s, s + epsilon_s);
    y_i ~ uniform(i - epsilon_i, i + epsilon_i);
    y_r ~ uniform(r - epsilon_r, r + epsilon_r);
  }
}
