function prepare_init()
    nc = netcdf('data/init.nc', 'c');
    nc{'beta'} = ncdouble();
    nc{'nu'} = ncdouble();
    nc{'sigma1'} = ncdouble();
    nc{'sigma2'} = ncdouble();
    nc{'beta'}(:) = 0.0033405;
    nc{'nu'}(:) = 0.61267;
    nc{'sigma1'}(:) = 0.0067957;
    nc{'sigma2'}(:) = 0.29879;
    ncclose(nc);
end
