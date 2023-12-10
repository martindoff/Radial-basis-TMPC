function g = RBF_cvx(x, beta, f_RBF_)
%RBF Compute g in f= g - h, where g, h are convex and where f is a radial 
% basis function approximation f = sum_i beta_i f_RBF_i(x) where f_RBF_i are
% RBF and beta_i coefficients. 
% Input:  - x: variable
%         - beta: coefficients
%         - f_RBF_: radial basis functions
% Output: - g: function g in f= g - h, where g, h are convex

N = length(f_RBF_); % number of RBF functions
g = 0;  % by convention, g includes the constant term 
for i=1:N
    if beta(i) >=0  % select only the RBF with positive coefficients
        g = g + beta(i)*f_RBF_{i}(x);
    end 
    
end 

end

