function f = RBF(x, beta, f_RBF_)
%RBF Radial basis function approximation f = sum_i beta_i f_RBF_i(x) where f_RBF_i are
% RBF and beta_i coefficients. 
% Input:  - x: variable
%         - beta: coefficients
%         - f_RBF_: radial basis functions
% Output: - f: approximation in terms of RBF

N = length(f_RBF_); % number of RBF functions
f = 0; 
for i=1:N
    f = f + beta(i)*f_RBF_{i}(x);
end 

end

