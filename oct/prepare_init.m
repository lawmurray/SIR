function prepare_init()
    nc = netcdf('data/init.nc', 'c');

    nc{'theta1_beta'} = ncdouble();
    nc{'theta2_beta'} = ncdouble();
    nc{'theta3_beta'} = ncdouble();
    nc{'theta1_nu'} = ncdouble();
    nc{'theta2_nu'} = ncdouble();
    nc{'theta3_nu'} = ncdouble();

    nc{'theta1_beta'}(:) = -6.42;
    nc{'theta2_beta'}(:) = 1.06;
    nc{'theta3_beta'}(:) = -1.42;
    nc{'theta1_nu'}(:) = -2.62;
    nc{'theta2_nu'}(:) = 3.33;
    nc{'theta3_nu'}(:) = 3.10;
    
    ncclose(nc);
end
