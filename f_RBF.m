function f_out = f_RBF(x, u, theta, c, rho, dtheta)
%F_RBF Compute RBF function approximation of f using sums of multiquadrics 
% Input:  - x: state 
%         - u: input
%         - theta: weights of RBF functions
%         - c: centers of RBFs
%         - rho: scalings of RBFs
%         - dtheta: parameter uncertainty rescaling
% Output: - f_out: RBF approximation

if nargin < 6
    dtheta = 1; 
end 
% Problem dimensions
[nx, N] = size(x);
%f_out = zeros(nx, N);

% Redefine input
%z = zeros(nx, N, nx);
z(:, :, 1) = [x(1, :); u];
z(:, :, 2) = x;

% Loop through state
for k=1:nx
    % Extract variables 
    theta_ = theta{k}*dtheta;
    c_ = c{k};
    rho_ = rho{k};
    N_RBF = length(theta_); 

    % Prepare input
    x_ = repmat(z(:,:,k), [1, 1, N_RBF]); 
    x_ = permute(x_, [1, 3, 2]); 

    % Multiquadric
    dx = repmat(sqrt(rho_)', [nx, 1, N]).*(x_ - repmat(c_, [1, 1, N]));
    dx = [ dx ; ones(1, N_RBF,N)];
    mult = norms(dx, 2, 1);


%     % Multiquadric
%     dx = x_ - repmat(c_, [1, 1, N]);
%     dx2 = sum(dx.^2, 1);
%     rdx = repmat(rho_', [1, 1, N]).*dx2;
%     mult = sqrt(1 + rdx);

    % RBF sum
    f = sum(repmat(theta_', [1, 1, N]).*mult, 2);
    
    % Reshape output
    f_out(k, :) = reshape(f,  1, N); 
end 
end

