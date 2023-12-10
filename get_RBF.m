function [f, g, h, theta, theta_g, theta_h, c_g, c_h, rho_g, rho_h, MAE_train] = get_RBF(N_samples, c_RBF, rho_RBF, p, input_train, y_train)
% GET_RBF Compute coefficients of the RBF approximation of a function and
% the coefficients of the DC decomposition
%

% Initialise
N_RBF = size(c_RBF,2); % number of RBF
f_RBF_ = cell(N_RBF, 1);
f = cell(p.nx, 1);
g = cell(p.nx, 1);
h = cell(p.nx, 1);
theta = cell(p.nx, 1);
theta_g = cell(p.nx, 1);
theta_h = cell(p.nx, 1);
c_g = cell(p.nx, 1);
c_h = cell(p.nx, 1);
rho_g = cell(p.nx, 1);
rho_h = cell(p.nx, 1);

% RBF function
for i=1:N_RBF
    f_RBF_{i} = @(x) multiquad(x, c_RBF(:, i), rho_RBF(i));  % RBF
end 

% Loop through each state 
for k=1:p.nx
    xi = input_train{k}; % input to RBF
    Y_train = y_train(k,:)';  % target
    X_train = zeros(N_samples,N_RBF);
    
    % Prepare training data
    for i=1:N_RBF
        for j=1:N_samples
            X_train(j, i) = f_RBF_{i}(xi(:, j));
        end 
    end
    
    % Least square
    theta{k} = X_train\Y_train;
    lambda = 10; 
    theta{k} = (X_train'*X_train + N_samples*lambda*eye(N_RBF))\X_train'*Y_train;

    % Decomposition
    theta_g{k} = theta{k}(theta{k} >= 0);
    c_g{k} = c_RBF(:, theta{k} >= 0);
    rho_g{k} = rho_RBF(theta{k} >= 0);

    theta_h{k} = -theta{k}(theta{k} < 0);
    c_h{k} = c_RBF(:, theta{k} < 0);
    rho_h{k} = rho_RBF(theta{k} < 0);

    
    % RBF approx
    f{k} = @(x) (RBF(x, theta{k}, f_RBF_));  % RBF approximation
    g{k} = @(x) (RBF_cvx(x, theta{k}, f_RBF_));  % g convex in f = g - h
    h{k} = @(x) (RBF_ccv(x, theta{k}, f_RBF_));  % h convex in f = g - h
    y_pred_train = f{k}(xi);  % prediction on training data
    
    % Fit evaluation 
    fprintf('Train evaluation for state %d \n', k)
    MAE_train(k, 1) =  mean(abs(y_pred_train - Y_train'))

end 
end 