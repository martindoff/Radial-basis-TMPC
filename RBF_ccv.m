function h = RBF_ccv(x, beta, f_RBF_)
%RBF Compute h in f= g - h, where g, h are convex and where f is a radial 
% basis function approximation f = sum_i beta_i f_RBF_i(x) where f_RBF_i are
% RBF and beta_i coefficients. 
% Input:  - x: variable
%         - beta: coefficients
%         - f_RBF_: radial basis functions
% Output: - h: function h in f= g - h, where g, h are convex

N = length(f_RBF_); % number of RBF functions
h = 0;  % by convention, h do not include the constant term beta(1)
for i=1:N
    if beta(i) < 0  % select only the RBF with negative coefficients
        h = h - beta(i)*f_RBF_{i}(x);
    end 
    
end 

end

