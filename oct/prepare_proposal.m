function prepare_proposal ()
    input_file = 'data/input.nc';
    output_file = 'results/posterior_bridge.nc';

    theta(1,:,1) = ncread(output_file, 'theta1_beta');
    theta(1,:,2) = ncread(output_file, 'theta2_beta');
    theta(1,:,3) = ncread(output_file, 'theta3_beta');
    theta(2,:,1) = ncread(output_file, 'theta1_nu');
    theta(2,:,2) = ncread(output_file, 'theta2_nu');
    theta(2,:,3) = ncread(output_file, 'theta3_nu');
    C(1,:,:) = cov(squeeze(theta(1,:,:)));
    C(2,:,:) = cov(squeeze(theta(2,:,:)));
    
    % write file
    try
        nccreate(input_file, 'C', 'Dimensions', {'f'; 2; 'n'; 3; 'n'; 3});
    catch
        % assume variables already exist...
    end
    ncwrite(input_file, 'C', C);
end
