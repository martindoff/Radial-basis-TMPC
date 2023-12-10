function [MAE_test] = test_fit(func, sqt_N_test, f, g, h, p, plt)
%TEST_FIT Test goodness of fit (works only for 2D problems, needs to be adapted to problem)

% generate test data
[X1_test, Y1_test] = meshgrid(linspace(p.x_min(1), p.x_max(1), sqt_N_test), linspace(p.u_min, p.u_max, sqt_N_test));
[X2_test, Y2_test] = meshgrid(linspace(p.x_min(1), p.x_max(1), sqt_N_test), linspace(p.x_min(2), p.x_max(2), sqt_N_test));
y_test = func([X2_test(:)'; Y2_test(:)'], Y1_test(:)', p);
input_test = {[X1_test(:)'; Y1_test(:)'], [X2_test(:)'; Y2_test(:)']}; % input to the RBF for each state
X_ = {X1_test, X2_test};
Y_ = {Y1_test, Y2_test};


% Loop through each state 
max_err = zeros(p.nx, 1);
MAE_test = zeros(p.nx, 1);
for k=1:p.nx
    xi_test = input_test{k};
    Y_test = y_test(k,:)';
    y_pred_test = f{k}(xi_test);  % prediction on test data
    
    % Fit evaluation 
    fprintf('Test evaluation for state %d \n', k)
    MAE_test(k) =  mean(abs(y_pred_test - Y_test'))
    %max_err(k) = max(abs(y_pred_test - Y_test'))


    % Plot results 
    if plt== true

        font_size = 15;
        line_size = 15;
        line_width = 2;
        
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
        view([-37.5 30])

    end 
    
end 

end

