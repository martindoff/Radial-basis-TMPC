function x = eul(f, x0, u, p, w)
% EUL Integrate dynamics using forward Euler method:
%     x[k+1] = x[k] + delta*(f(x[k], u[k]))
% Input:  - f: function to integrate
%         - x0: initial condition
%         - u: control input
%         - p: structure of parameters
%         - w: external disturbance
% Output: - x: trajectory

if nargin < 5
  w = 0; % default zero disturbance (nominal case) 
end

N = size(u, 2);
N_state = size(x0, 1); 
x = zeros(N_state, N+1);
x(:, 1) = x0;
for i=1:N
    x(:, i+1) = x(:, i) + p.delta*(f(x(:, i), u(:, i), p) + w);
end 
end

