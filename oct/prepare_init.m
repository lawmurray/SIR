function prepare_init()
    nc = netcdf('data/init.nc', 'c');
    nc{'beta'} = ncdouble();
    nc{'nu'} = ncdouble();
    nc{'sigma1'} = ncdouble();
    nc{'sigma2'} = ncdouble();
    nc{'beta'}(:) = 1.663/763;
    nc{'nu'}(:) = 0.441;
    nc{'sigma1'}(:) = 100.0;
    nc{'sigma2'}(:) = 100.0;
    ncclose(nc);
end
