clear all; close all; clc

%% Simulation parameters
p = param_init();                  % initialise problem parameters 
T_sim = 35;                        % simulation time
Tmax = 10;                         % MPC horizon
p.delta = T_sim/Tmax;              % time step
maxiter = 1;
Q = [0 0;0 .05]; R = .01;          % cost matrices
sqrtQ = Q^0.5; sqrtR = R^0.5; 
dtheta = 1; p.dtheta = dtheta;

%% Initialise problem
x = zeros(p.nx, Tmax+1);           % state
u = zeros(p.nu, Tmax);             % input
x(:,1) = p.x_init;                 % initial state
u_0 = p.u_init*ones(p.nu, Tmax);   % initial input
x_r = p.h_r.*ones(p.nx, Tmax+1);   % reference state
u_r = p.u_r.*ones(p.nu, Tmax);     % reference input
t_avg= 0;                          % counter for average time / iteration
iter_count = 0;                    % counter for number of iterations

%% DC decomposition
% Define RBF 
N_samples = 1000; 
sqt_N_RBF = 3;
N_RBF = sqt_N_RBF^2;               % number of RBF (only works in 2D)
[X1_RBF, X2_RBF] = meshgrid(linspace(p.x_min(1)-3, p.x_max(1)+3, sqt_N_RBF),...
                            linspace(p.x_min(1)-3, p.x_max(1)+3, sqt_N_RBF));
c_RBF = [X1_RBF(:)';X2_RBF(:)'];   % RBF centers
rho_RBF = ones(N_RBF);             % RBF scalings 

% Generate training data
[x_train, u_train, y_train] = gen_train(@dynamics, N_samples, p);
input_train = {[x_train(1, :); u_train], x_train}; 
% input to the RBF for each state (needs adaptation to specific dynamics)
  
% Fit RBF and get decomposition
[f_RBF_, g_RBF, h_RBF, theta,... 
 theta_g, theta_h, c_g, c_h,...
 rho_g, rho_h, MAE_train] = get_RBF(N_samples, c_RBF, rho_RBF, ...
                                    p, input_train, y_train); 

% Function wrappers (problem specific)
f = @(x, u, p) ([f_RBF_{1}([x(1, :); u]); f_RBF_{2}(x)]);
g_ = @(x, u, p) ([g_RBF{1}([x(1, :); u]); g_RBF{2}(x)]);
h_ = @(x, u, p) ([h_RBF{1}([x(1, :); u]); h_RBF{2}(x)]);

% Test fit (will only work for specific coupled tank problem)
plt = true; 
dim_N_test = 10; % test points per dimension
MAE = test_fit(@dynamics, dim_N_test, f_RBF_, g_RBF, h_RBF, p, plt); 
p.w_max = MAE_train;

%% Terminal set
% Linearise g at reference
[A1_term, B1_term] = linearise(p.h_r, p.u_r, theta_g, c_g, rho_g, p);

 % Linearise h at reference 
[A2_term, B2_term] = linearise(p.h_r, p.u_r, theta_h, c_h, rho_h, p); 

% Discretise
A_d_term = eye(p.nx) + p.delta*(A1_term - A2_term);
B_d_term = p.delta*(B1_term - B2_term);

% Set optimisation
[K,P,V,gam,beta] = term_comp(A_d_term, B_d_term,Q,R,p)
sqrtP = P^0.5; sqrtV = V^0.5; 

%% Feasible trajectory
x_0 = zeros(p.nx, Tmax+1);
x_0(:,1) = p.x_init;
for i=1:Tmax
    % Control input (dummy controller)
    u_0(:,i) = max(min(K*(x_0(:,i)-x_r(:,i)) + u_r(:,i), p.u_max), p.u_min);
    
    % Generate feasible trajectory
    x_0(:,i+1) = x_0(:,i) + p.delta*(f_RBF(x_0(:,i), u_0(:,i), ...
     theta_g, c_g, rho_g) - f_RBF(x_0(:,i), u_0(:,i), theta_h, c_h, rho_h));
end 

% initial seed for the nominal input perturbation
c_0 = u_0 - K*x_0(:, 1:end-1);

%% Linearise system (CT)
[A1, B1] = linearise(x_0, u_0, theta_g, c_g, rho_g, p);  % linearise g
[A2, B2] = linearise(x_0, u_0, theta_h, c_h, rho_h, p);  % linearise h

% Closed loop 
Phi1 = A1 + B1.*K;
Phi2 = A2 + B2.*K;

% Trajectory in phase plot
figure()
hold on 
plot(x_0(1, :), x_0(2, :))
rectangle('Position',[0 0 p.x_max(1) p.x_max(2)])
axis([p.x_min(1)-1 p.x_max(1)+1 p.x_min(2)-1 p.x_max(2)+1])
ellipse(V, x_r(:, end))
xlabel('$x_1$', 'Interpreter','latex')
ylabel('$x_2$', 'Interpreter','latex')
legend({'$x_0$', 'terminal set'}, 'Interpreter','latex')
grid on
                            
%% MPC loop
tol = 1e-8; delta_c = tol*10;
for t = 1:Tmax
    fprintf("************** Solving problem at time %d/%d ******************\n", t, Tmax)

    iter = 1; c_old = 0; 
    while (iter <= maxiter && norm(delta_c) >= tol)
        fprintf("iteration %d/%d (time=%d)\n", iter, maxiter, t)
       
        %% Optimisation
        tic 
        [c, x_lb, x_ub, cvx_optval] = cvx_optimisation(x_0, u_0, c_0, ...
                                      x_r, u_r, Phi1, Phi2, B1, B2, K, ...
                                      sqrtR, sqrtQ, sqrtV, sqrtP, p, ...
                                      theta_g, c_g, rho_g, theta_h, c_h, rho_h); 
         
        t_elasped = toc 
        t_avg = t_avg + t_elasped;  
        iter_count = iter_count+1; 
        
        %% Update feasible trajectory
        % Update of nominal input disturbance
        c_0 = c_0 + c; 

         % State update
         x_ = x_0(:, 1);
         for i=1:Tmax
             u_0(:,i) = c_0(:,i) + K*x_; 
             x_0(:,i+1) = x_0(:,i) + p.delta*(f_RBF(x_0(:,i), u_0(:,i), ...
                          theta_g, c_g, rho_g) - f_RBF(x_0(:,i), ...
                          u_0(:,i), theta_h, c_h, rho_h));
             x_ = x_0(:,i+1);
         end 
         
         % Linearise system (CT)
         [A1, B1] = linearise(x_0, u_0, theta_g, c_g, rho_g, p);  % linearise g
         [A2, B2] = linearise(x_0, u_0, theta_h, c_h, rho_h, p);  % linearise h

         % Closed loop
         Phi1 = A1 + B1.*K;
         Phi2 = A2 + B2.*K;

        % Cost and optimal solution
        J = cvx_optval; 
        delta_c = c-c_old;
        c_old = c;
        J_iter(iter) = J;

       
        iter = iter + 1;
    end  

    Iter(t) = iter;
    J_time(t) = J;

    %% Update system 
    % Control input
    u_opt(t) = K*x(:,t) + c_0(1);

    % State update
    w_true = rand*2*p.w_max-p.w_max;
    x(:,t+1) = x(:,t) + p.delta*dynamics(x(:,t), u_opt(t), p);
    %x(:,t+1) = x(:,t) + p.delta*(f_RBF(x(:,t), u_opt(t), theta_g, ...
    % c_g, rho_g) - f_RBF(x(:,t), u_opt(t), theta_h, c_h, rho_h) + w_true);
    
    % Nominal trajectory update
    c_0 = [c_0(2:Tmax) u_r(:, Tmax)-K*x_r(:, Tmax)];
    u_0(1:Tmax-1) = u_0(2:Tmax);
    x_0(:, 1) = x(:, t+1);
    for i=1:Tmax-1
        x_0(:,i+1) = x_0(:,i) + p.delta*(f_RBF(x_0(:,i), u_0(:,i), theta_g, ...
            c_g, rho_g) - f_RBF(x_0(:,i), u_0(:,i), theta_h, c_h, rho_h));
    end 

    % Dual mode paradigm
    u_0(:, Tmax) = K*(x_0(:, Tmax) - x_r(:, Tmax)) + u_r(:, Tmax); 
    x_0(:,Tmax+1) = x_0(:,Tmax) + p.delta*(f_RBF(x_0(:,Tmax), u_0(:,Tmax), ...
    theta_g, c_g, rho_g) - f_RBF(x_0(:,Tmax), u_0(:,Tmax), theta_h, c_h, rho_h));
    
    % Compute Linearization
    [A1, B1] = linearise(x_0, u_0, theta_g, c_g, rho_g, p);  % linearise g
    [A2, B2] = linearise(x_0, u_0, theta_h, c_h, rho_h, p);  % linearise h
         
    % Closed loop
    Phi1 = A1 + B1.*K;
    Phi2 = A2 + B2.*K;
    
end

fprintf('Average time per iteration of the optimisation: %.2f s\n', t_avg/iter_count)

%% Plot results
t = (0:Tmax)*p.delta;
figure
subplot(3,1,1)
stairs(t, [u_opt u_opt(end)])
axis([0 T_sim min(u_opt)/1.2 max(u_opt)*1.2])
ylabel('Control, $u_t$ (V)', 'Interpreter','latex')
grid on
subplot(3,1,2)
stairs(t,x(1,:))
ylabel('State, $[x_t]_1$ (cm)', 'Interpreter','latex')
grid on
subplot(3,1,3)
stairs(t,x(2,:))
hold on
plot(t,x_r(2,:),'--')
axis([0 T_sim 0 20])
ylabel('State, $[x_t]_2$ (cm)', 'Interpreter','latex')
xlabel('time step, t', 'Interpreter','latex')
grid on