function [c, xlb, xub, cvx_optval] = cvx_optimisation(x_0, u_0, c_0, x_r, u_r, Phi1, Phi2, B1, B2, K, sqrtR, sqrtQ, sqrtV, sqrtP, p, theta_g, c_g, rho_g, theta_h, c_h, rho_h)
%CVX_OPTIMISATION Solve optimisation problem

% Initialisation
[nu, N] = size(u_0);    % number of inputs / horizon 
nx = size(x_0, 1);      % number of states       
nv = 2^nx;              % number of vertices
x_r_term = x_r(:,end);
x_r = x_r(:, 1:end-1);
x_0 = x_0(:, 1:end-1);

% Constraint sets
U_min = p.u_min; 
U_max = p.u_max;
X_min = p.x_min.*ones(nx, N);
X_max = p.x_max.*ones(nx, N);

% Some useful matrices
Phi1_c =num2cell(Phi1,[1,2]); 
Phi1_ = blkdiag(Phi1_c{:});
Phi2_c =num2cell(Phi2,[1,2]); 
Phi2_ = blkdiag(Phi2_c{:});
B1_c =num2cell(B1,[1,2]);
B1_ = blkdiag(B1_c{:});
B2_c =num2cell(B2,[1,2]);
B2_ = blkdiag(B2_c{:});

cvx_begin quiet
  cvx_solver mosek
   variables l_x(N) l_u(N) l_n(1) xub(nx, N+1) xlb(nx, N+1) c(nu, N)
   expressions x_vertex(nx, N+1, nv) x(nx, N+1);
   minimize(sum(l_x + l_u) + l_n)
   subject to
   
   % Define vertices
   x_vertex(:, :, 1) = [xlb(1, :); xlb(2, :)];
   x_vertex(:, :, 2) = [xub(1, :); xlb(2, :)];
   x_vertex(:, :, 3) = [xlb(1, :); xub(2, :)];
   x_vertex(:, :, 4) = [xub(1, :); xub(2, :)];

   for l=1:nv
       % Current vertex 
       x = x_vertex(:, 1:end-1, l);
       x_term = x_vertex(:, end, l);

       % Useful variables
       dx = reshape(x-x_0, [nx*N, 1]);
       dc = reshape(c, [nu*N, 1]);
       Phi1_dx = reshape(Phi1_ * dx, [nx, N]);
       Phi2_dx = reshape(Phi2_ * dx, [nx, N]);
       B1_c = reshape(B1_ * dc, [nx, N]);
       B2_c = reshape(B2_ * dc, [nx, N]);

       % Objective
       norm(sqrtQ*(x - x_r), 1) <= l_x;
       norm(sqrtR*(K*x + c + c_0 - u_r), 1) <= l_u;
       norm(sqrtP*(x_term - x_r_term)) <= l_n;

       % Input constraints
       U_min <= K*x + c + c_0;
       U_max >= K*x + c + c_0;

       % State constraints
       X_min <= x;
       X_max >= x;

       % Initial conditions
       x(:, 1) == x_0(:, 1); 

       % Tube constraint
       xlb(:, 2:end) <= x + p.delta*(f_RBF(x_0, u_0, theta_g, c_g, rho_g, 1/p.dtheta) + Phi1_dx + B1_c ...
           -f_RBF(x, K*x + c + c_0, theta_h, c_h, rho_h, p.dtheta) - p.w_max.*ones(nx, N) )
       xub(:, 2:end) >= x + p.delta*(f_RBF(x, K*x + c + c_0, theta_g, c_g, rho_g, p.dtheta) ...
           - (f_RBF(x_0, u_0, theta_h, c_h, rho_h, 1/p.dtheta) + Phi2_dx + B2_c) + p.w_max.*ones(nx, N) )

       % Terminal set
       norm(sqrtV*(x_term - x_r_term)) <= 1;
   end 

cvx_end
end