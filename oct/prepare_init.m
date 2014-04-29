function prepare_init()
    init_file = 'data/init.nc';
    
    % create file
    try
        nccreate(init_file, 'theta_beta', 'Dimensions', {'n'; 3});
        nccreate(init_file, 'theta_nu', 'Dimensions', {'n'; 3});
    catch
        % assume variables already exist...
    end
    
    % write file
    ncwrite(init_file, 'theta_beta', [-6.42 1.06 -1.42]');
    ncwrite(init_file, 'theta_nu', [-2.62 3.33 3.10]');
end
