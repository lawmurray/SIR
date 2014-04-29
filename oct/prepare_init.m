function prepare_init()
    init_file = 'data/init.nc';
    
    % write file
    try
        nccreate(init_file, 'theta', 'Dimensions', {'f'; 2; 'n'; 3});
    catch
        % assume variable already exists...
    end
    ncwrite(init_file, 'theta', [-18.203 3.105 -0.150; -49.290 38.744 13.305]);
end
