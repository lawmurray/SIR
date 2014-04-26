function prepare_init()
    nc = netcdf('data/init.nc', 'c');

    nc{'theta1_beta'} = ncdouble();
    nc{'theta2_beta'} = ncdouble();
    nc{'theta3_beta'} = ncdouble();
    nc{'theta1_nu'} = ncdouble();
    nc{'theta2_nu'} = ncdouble();
    nc{'theta3_nu'} = ncdouble();

    nc{'theta1_beta'}(:) = log(0.003);
    nc{'theta2_beta'}(:) = 1.0;
    nc{'theta3_beta'}(:) = 0.5;
    nc{'theta1_nu'}(:) = log(0.6);
    nc{'theta2_nu'}(:) = 1.0;
    nc{'theta3_nu'}(:) = 0.5;
    
    ncclose(nc);
end
