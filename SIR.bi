/**
 * SIR model.
 */
model SIR {
  const h = 1.0e-3; // time step
  const eps = 1.0e-1;
  
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
  const lambda = 2; // tempering

  sub parameter {
    beta ~ uniform(0.0, 1.0);
    nu ~ uniform(0.0, 1.0);
    sigma1 ~ uniform(0.0, 1000.0);
    sigma2 ~ uniform(0.0, 1000.0);
  }

  sub proposal_parameter {
    //beta ~ truncated_gaussian(beta, 1.0e-5, 0.0, 1.0);
    //nu ~ truncated_gaussian(nu, 1.0e-3, 0.0, 1.0);
    //sigma1 ~ truncated_gaussian(sigma1, 1.0, 0.0, 1000.0);
    //sigma2 ~ truncated_gaussian(sigma2, 1.0, 0.0, 1000.0);
  }

  sub initial {
    s <- s0;
    i <- i0;
    r <- r0;
  }

  sub transition(delta = h) {
    w1 ~ wiener();
    w2 ~ wiener();

    ode(h = h, alg = 'RK4') {
      ds/dt = -beta*i*s + sigma1*s*w1;
      di/dt = beta*i*s - nu*i - sigma1*s*w1 + sigma2*i*w2;
      dr/dt = nu*i - sigma2*i*w2;
    }
  }

  sub bridge {
    inline s_k = s_sf2*exp(-0.5*(t_next_obs - t_now)**2/s_ell2);
    inline s_mu = log(s)*s_k/s_sf2;
    inline s_sigma = sqrt(s_sf2 - s_k*s_k/s_sf2 + eps*eps);

    y_s ~ log_normal(s_mu, s_sigma);

    inline i_k = i_sf2*exp(-0.5*(t_next_obs - t_now)**2/i_ell2);
    inline i_mu = log(i)*i_k/i_sf2;
    inline i_sigma = sqrt(i_sf2 - i_k*i_k/i_sf2 + eps*eps);

    y_i ~ log_normal(i_mu, i_sigma);

    inline r_k = r_sf2*exp(-0.5*(t_next_obs - t_now)**2/r_ell2);
    inline r_mu = log(r)*r_k/r_sf2;
    inline r_sigma = sqrt(r_sf2 - r_k*r_k/r_sf2 + eps*eps);

    y_r ~ log_normal(r_mu, r_sigma);
  }

  sub observation {
    y_s ~ gaussian(s, eps);
    y_i ~ gaussian(i, eps);
    y_r ~ gaussian(r, eps);
  }
}
