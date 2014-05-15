function prepare_init()
    init_file = 'data/init.nc';
    
    % write file
    try
        nccreate(init_file, 'theta', 'Dimensions', {'f'; 2; 'n'; 3});
    catch
        % assume variable already exists...
    end
    ncwrite(init_file, 'theta', [-19.203 3.049 1.044; -1.029 2.305 0.989]);
end
