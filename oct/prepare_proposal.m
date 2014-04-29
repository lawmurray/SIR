function prepare_proposal ()
    input_file = 'data/input.nc';
    output_file = 'results/posterior_bridge.nc';

    theta_beta(:,1) = ncread(output_file, 'theta1_beta');
    theta_beta(:,2) = ncread(output_file, 'theta2_beta');
    theta_beta(:,3) = ncread(output_file, 'theta3_beta');
    U_beta = chol(cov(theta_beta));
    
    theta_nu(:,1) = ncread(output_file, 'theta1_nu');
    theta_nu(:,2) = ncread(output_file, 'theta2_nu');
    theta_nu(:,3) = ncread(output_file, 'theta3_nu');
    U_nu = chol(cov(theta_nu));
    
    % create file
    try
        nccreate(input_file, 'U_beta', 'Dimensions', {'n'; 3; 'n'; 3});
        nccreate(input_file, 'U_nu', 'Dimensions', {'n'; 3; 'n'; 3});
    catch
        % assume variables already exist...
    end
    
    % write file
    ncwrite(input_file, 'U_beta', U_beta);
    ncwrite(input_file, 'U_nu', U_nu);
    
    theta(:,4) = ncread(file, 'theta1_nu');
    theta(:,5) = ncread(file, 'theta2_nu');
    theta(:,6) = ncread(file, 'theta3_nu');
end
