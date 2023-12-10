function [outputArg1,outputArg2] = ellipse(P, x)
%ELLIPSE plots an ellipse of the form xPx = 1

if nargin < 2
     x = [0;0]; 
end 

R = chol(P);
t = linspace(0, 2*pi, 100); % or any high number to make curve smooth
z = [cos(t); sin(t)];
e = inv(R) * z;
plot(e(1,:)+x(1), e(2,:)+x(2),'-r', 'LineWidth', 2)

end

