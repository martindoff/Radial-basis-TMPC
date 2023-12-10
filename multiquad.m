function y = multiquad(x, c, rho)
%MULTIQUAD multiquadric radial basis function
% Input: - x: input
%        - c: center
%        - rho: scaling

dx = x - c;
y = sqrt(1 + rho*dot(dx, dx)); 
end

