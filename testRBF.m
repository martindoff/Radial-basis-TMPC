%% RBF test
clear all; close all; clc
% Initialise problem
p = param_init();

% generate training data
N_samples = 1000;
x_train = (p.x_max-p.x_min).*rand(p.nx, N_samples) + p.x_min;
u_train = (p.u_max-p.u_min).*rand(p.nu, N_samples) + p.u_min;
y_train = dynamics(x_train, u_train, p);
input_train = {[x_train(1, :); u_train], x_train}; % input to the RBF for each state

% generate test data
sqt_N_test = 10;
N_test = sqt_N_test^2;
[X1_test, Y1_test] = meshgrid(linspace(p.x_min(1), p.x_max(1), sqt_N_test), linspace(p.u_min, p.u_max, sqt_N_test));
[X2_test, Y2_test] = meshgrid(linspace(p.x_min(1), p.x_max(1), sqt_N_test), linspace(p.x_min(2), p.x_max(2), sqt_N_test));
%x_test = [linspace(p.x_min(1), p.x_max(1), N_test); linspace(p.x_min(2), p.x_max(2), N_test)]; 
%u_test = linspace(p.u_min, p.u_max, N_test);
y_test = dynamics([X2_test(:)'; Y2_test(:)'], Y1_test(:)', p);
input_test = {[X1_test(:)'; Y1_test(:)'], [X2_test(:)'; Y2_test(:)']}; % input to the RBF for each state
X_ = {X1_test, X2_test};
Y_ = {Y1_test, Y2_test};
% Normalise training data
%[x_train, avg_x, std_x] = scale(x_train);


% Load the training data.
%[x,t] = simplefit_dataset;
%N_samples = length(x);

% Define RBF 
sqt_N_RBF = 4;  % number of RBF per dimension 
N_RBF = sqt_N_RBF^2; % number of RBF
[X1_RBF, X2_RBF] = meshgrid(linspace(p.x_min(1), p.x_max(1), sqt_N_RBF), linspace(p.x_min(1), p.x_max(1), sqt_N_RBF)); % RBF centers
c_RBF = [X1_RBF(:)';X2_RBF(:)'];

rho_RBF = ones(size(c_RBF)); % RBF scalings 
f_RBF = cell(N_RBF, 1);

for i=1:N_RBF
    f_RBF{i} = @(x) multiquad(x, c_RBF(:, i), rho_RBF(i));  % RBF
end 

% Loop through each state 
f = cell(p.nx, 1);
for k=1:p.nx
    xi = input_train{k}; % input to RBF
    xi_test = input_test{k};
    Y_train = y_train(k,:)';  % target
    Y_test = y_test(k,:)';
    X_train = zeros(N_samples,N_RBF);
    X_test = zeros(N_test,N_RBF);
    
    % Prepare training data
    for i=1:N_RBF
        for j=1:N_samples
            X_train(j, i) = f_RBF{i}(xi(:, j));
        end 
    end

    % Prepare test data
    for i=1:N_RBF
        for j=1:N_test
            X_test(j, i) = f_RBF{i}(xi_test(:, j));
        end 
    end
    
    % Least square
    beta = X_train\Y_train;
    %beta = beta + randn(size(beta)).*beta*0.1
    %lambda = 0;
    %I = eye(N_RBF+1);
    %beta = (X_train'*X_train + N_samples*lambda*I)\(X_train'*Y_train);

    % 
    % cvx_begin % Start CVX 
    %  cvx_precision high
    %  cvx_solver mosek
    %     variables beta(N_RBF+1)
    %     minimize((X_train*beta - Y_train)'*(X_train*beta - Y_train)) %last term is for windmilling (maximise final battery SOC)
    %     subject to
    % cvx_end

    
    % RBF approx
    f{k} = @(x) (RBF(x, beta, f_RBF));  % RBF approximation
    g{k} = @(x) (RBF_cvx(x, beta, f_RBF));  % g convex in f = g - h
    h{k} = @(x) (RBF_ccv(x, beta, f_RBF));  % g convex in f = g - h
    y_pred_train = f{k}(xi);  % prediction on training data
    y_pred_test = f{k}(xi_test);  % prediction on test data
    % Fit evaluation 

    fprintf('Evaluation for state %d \n', k)
    MAE_train =  mean(abs(y_pred_train - Y_train'))
    MAE_test =  mean(abs(y_pred_test - Y_test'))

    % Plot results 
    font_size = 15;
    line_size = 15;
    line_width = 2;
    
    % Scatter plot
%     figure
%     hold on
%     scatter3(xi_test(1,:), xi_test(2,:), Y_test, '+r')
%     scatter3(xi_test(1,:), xi_test(2,:), g{k}(xi_test)-h{k}(xi_test), '.b')
%     scatter3(xi_test(1,:), xi_test(2,:), g{k}(xi_test), '.g')
%     scatter3(xi_test(1,:), xi_test(2,:), h{k}(xi_test), '.black')
%     legend('data', 'g-h RBF', 'g RBF', 'h RBF')


    % Create grid data
    X = X_{k};
    Y = Y_{k};
    F = zeros(size(X));
    G = zeros(size(X));
    H = zeros(size(X));

    for i=1:sqt_N_test
        for j=1:sqt_N_test
            in = [X(i, j); Y(i, j)];
            F(i, j) = g{k}(in)-h{k}(in);
            G(i, j) = g{k}(in);
            H(i, j) = h{k}(in);
        end 
    end 
    
    % Surface plot
    figure
    hold on
    scatter3(xi_test(1,:), xi_test(2,:), Y_test, '+r','Linewidth',line_width)
    surf(X, Y, F, 'FaceAlpha', 0.5,'FaceColor',[0 0 1])
    surf(X, Y, G, 'FaceAlpha', 0.5,'FaceColor',[1 0 0])
    surf(X, Y, H, 'FaceAlpha', 0.5,'FaceColor',[0 1 0])
    legend('data', 'RBF: g-h', 'RBF: g', 'RBF: h', 'fontsize',font_size,'Interpreter','latex')
    xlabel('$x_1$','fontsize',font_size,'Interpreter','latex')
    ylabel('$x_2$','fontsize',font_size,'Interpreter','latex')
    zlabel('$f(x_1, x_2)$','fontsize',font_size,'Interpreter','latex')
    set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
    grid on
end 