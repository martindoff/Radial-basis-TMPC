function [x_train, u_train, y_train] = gen_train(func, N_samples, p)
%GEN_TRAIN Generate N_samples training data points
% Input:  - func: function to generate training output
%         - N_samples: number of training points
%         - p: structure of parameters
% Output: - x_train, u_train: training input
%         - y_train: training output
x_train = (p.x_max-p.x_min).*rand(p.nx, N_samples) + p.x_min;
u_train = (p.u_max-p.u_min).*rand(p.nu, N_samples) + p.u_min;
y_train = func(x_train, u_train, p);

% Normalise training data
%[x_train, avg_x, std_x] = scale(x_train);

end

