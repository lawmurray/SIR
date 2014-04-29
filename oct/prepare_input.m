function prepare_input ()
  obs_file = 'data/obs.nc';
  input_file = 'data/input.nc';
    
  % data
  t = ncread(obs_file, 'time');
  y = log(ncread(obs_file, 'y_i'));

  % fit
  inffunc = @infExact;
  meanfunc = @meanZero;
  hyp.mean = [];
  covfunc = @covSEiso;
  hyp.cov = log([5.0; 5.0]);
  likfunc = @likGauss;
  hyp.lik = log(1.0e-2);

  hyp = minimize(hyp, @gpwrap, -1000, @infExact, meanfunc, covfunc, likfunc, t, y);

  ell2 = exp(2*hyp.cov(1));
  sf2 = exp(2*hyp.cov(2));

  % write weight function parameters to input file
  try
      nccreate(input_file, 'i_ell2');
      nccreate(input_file, 'i_sf2');
      nccreate(input_file, 's0');
      nccreate(input_file, 'i0');
      nccreate(input_file, 'r0');
  catch
      % assume variables already exist
  end
  ncwrite(input_file, 'i_ell2', ell2);
  ncwrite(input_file, 'i_sf2', sf2);
  
  % write initial conditions to input file
  ncwrite(input_file, 's0', 760);
  ncwrite(input_file, 'i0', 3);
  ncwrite(input_file, 'r0', 0);
end
