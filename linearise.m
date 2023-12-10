function [A, B] = linearise(x_0, u_0, theta, c, rho, p)
% LINEARISE Compute continuous-time linearised model of a RBF model around 
% a guess trajectory. Noting f(x, u) = theta'f_RBF(x, u), the function 
% computes the A, B continuous-time linear model around (x_0, u_0)
%
%
% Input:  - x_0: guess state trajectory
%         - u_0: guess input trajectory
%         - theta: coefficients of the RBF
%         - c: vector of centers of the RBF 
%         - rho: vector of scalings of the RBF
%         - p: structure of parameters 
% Output: - A, B: continuous-time linear model of f
 

% Problem dimensions 
[N_input, N] = size(u_0);
[N_state, ~] = size(x_0);

% Cap last state
if N > 1
    x_0 = x_0(:, 1:end-1);
end

% Linearised continuous-time model
A = zeros(N_state, N_state, N);
B = zeros(N_state, N_input, N);

a11=0; b11=0; a21=0; a22=0;  

% x1
N_RBF = length(theta{1}); 
for i=1:N_RBF
    xi_0 = x_0(1, :) - c{1}(1, i);
    ui_0 = u_0 - c{1}(2, i);
    ri_0 = xi_0.^2 + ui_0.^2;

   
    a11 = a11 + xi_0 * rho{1}(i) * theta{1}(i)./sqrt(1+rho{1}(i)*ri_0);
    b11 = b11 + ui_0 * rho{1}(i) * theta{1}(i)./sqrt(1+rho{1}(i)*ri_0);
end

% x2
N_RBF = length(theta{2}); 
for i=1:N_RBF
    xi_0 = x_0 - c{2}(:, i); 
    ri_0 = xi_0(1, :).^2 + xi_0(2, :).^2; 
    a21 = a21 + xi_0(1, :) * rho{2}(i) * theta{2}(i)./sqrt(1+rho{2}(i)*ri_0);
    a22 = a22 + xi_0(2, :) * rho{2}(i) * theta{2}(i)./sqrt(1+rho{2}(i)*ri_0);
end

A(1, 1, :) = a11; A(2, 1, :) = a21; A(2, 2, :) = a22; B(1, 1, :) = b11;
 
end

