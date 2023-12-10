function p = param_init()
%PARAM_INIT Summary of this function goes here
%   Detailed explanation goes here

p.g = 981;                              % gravity acceleration (cm/s^2)
p.k_p = 3.3;                            % pump gain (cm3 s-1 V-1)
d = 4.4;                                % tank diameter (cm3)
p.A = d^2*pi/4;                         % tank area (cm2)
p.A1 = sqrt((p.k_p*7.3)^2/(2*p.g*16));  % sigma_1*a_1 (obtained experimentally)
p.A2 = sqrt((p.k_p*7.3)^2/(2*p.g*15));  % sigma_2*a_2 (obtained experimentally)
h2_r = 15;                              % reference height tank 2
h1_r = (p.A2/p.A1)^2*h2_r;              % reference height tank 1
p.h_r = [h1_r; h2_r];                   % reference state
p.u_r =p.A1/p.k_p*sqrt(2*p.g*h1_r);     % reference input (6.1 - 9.3)
p.u_init = p.u_r;                       % initial input
p.x_init = [2; 1];                      % initial state
p.u_min = 0;                            % min input
p.u_max = 24;                           % max input
p.x_min = [0.1; 0.1];                   % min state
p.x_max = [30; 30];                     % max state
p.u_term = 5;                           % terminal set bound on input 
p.x_term = 5;                           % terminal set bound on state
p.nx = length(p.x_min);                 % number of states
p.nu = length(p.u_min);                 % number of inputs
p.delta = 0.3;                          % time step
end

