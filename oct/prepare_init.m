function prepare_init()
    init_file = 'data/init.nc';
    
    % write file
    try
        nccreate(init_file, 'theta', 'Dimensions', {'f'; 2; 'n'; 3});
    catch
        % assume variable already exists...
    end
    ncwrite(init_file, 'theta', [-6.42 1.06 -1.42; -2.62 3.33 3.10]);
end
