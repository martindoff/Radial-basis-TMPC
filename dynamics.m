function x_ = dynamics(x, u, param)
%DYNAMICS compute the dynamic update x_ = f(x, u)
%   Input:  - x: state at step k
%           - u: input at step k
%           - param: structure of parameters
%   Output: - x_: state at step k+1

f1 = (param.k_p*u - param.A1*sqrt(2*param.g*x(1, :)))/param.A;
f2 = (param.A1*sqrt(2*param.g*x(1, :)) - param.A2*sqrt(2*param.g*x(2, :)))/param.A;

x_ = [f1; f2]; %x_ = x + param.delta*[f1; f2];
end

