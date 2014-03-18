% Copyright (C) 2013
% Author: Lawrence Murray <lawrence.murray@csiro.au>
function prepare_input ()
  % data
  nc = netcdf('data/obs.nc', 'r');
  t = nc{'time'}(:);
  y = log(nc{'y_i'}(:));
  ncclose(nc);

  % fit
  inffunc = @infExact;
  meanfunc = @meanZero;
  hyp.mean = [];
  covfunc = @covSEiso;
  hyp.cov = log([5.0; 1.0]);
  likfunc = @likGauss;
  hyp.lik = log(1.0e-3);

  hyp = minimize(hyp, @gpwrap, -1000, @infExact, meanfunc, covfunc, likfunc, t, y);

  ell2 = exp(2*hyp.cov(1));
  sf2 = exp(2*hyp.cov(2));

  % write weight function parameters to input file
  nc = netcdf('data/input.nc', 'w');
  %nc{'i_ell2'} = ncdouble();
  %nc{'i_sf2'} = ncdouble();
  nc{'i_ell2'}(:) = ell2;
  nc{'i_sf2'}(:) = sf2;
  ncclose(nc);
end
