function prepare_proposal ()
    init_file = 'data/init.nc';
    input_file = 'data/input.nc';
    output_file = 'results/posterior.nc';

    theta(1,:,:) = squeeze(ncread(output_file, 'theta')(:,1,:));
    theta(2,:,:) = squeeze(ncread(output_file, 'theta')(:,2,:));
    C(1,:,:) = cov(squeeze(theta(1,:,:)));
    C(2,:,:) = cov(squeeze(theta(2,:,:)));
    
    % write input file
    try
        nccreate(input_file, 'C', 'Dimensions', {'f'; 2; 'n'; 3; 'n'; 3});
    catch
        % assume variables already exist...
    end
    ncwrite(input_file, 'C', C);
    
    % write init file
    theta = squeeze(ncread(output_file, 'theta')(end,:,:));
    try
        nccreate(init_file, 'theta', 'Dimensions', {'f'; 2; 'n'; 3});
    catch
        % assume variables already exist...
    end
    ncwrite(init_file, 'theta', theta);
end
