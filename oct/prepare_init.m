function prepare_init()
    init_file = 'data/init.nc';
    
    % write file
    try
        nccreate(init_file, 'theta', 'Dimensions', {'f'; 2; 'n'; 3});
    catch
        % assume variable already exists...
    end
    ncwrite(init_file, 'theta', [-6.20 1.00 1.13; -1.64 1.00 4.32]);
end
