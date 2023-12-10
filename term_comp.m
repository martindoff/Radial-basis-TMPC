function [K,P,V, gamma, bet] = term_comp(A,B,Q,R,p)

A_term =A; B_term = B; pp.Q = Q; pp.R = R; C = Q(end, :);
G=[1/p.u_term; -1/p.u_term]; F=zeros(2,2); h=[1;1];

n = size(A_term,2);
m = size(B_term,2);

w(1:2,1) = [p.w_max(1);p.w_max(2)];
w(1:2,2) = [-p.w_max(1);-p.w_max(2)];
w(1:2,3) = [p.w_max(1);-p.w_max(2)];
w(1:2,4) = [-p.w_max(1);p.w_max(2)];
w = w*p.delta;

% Maximize invariant feasible ellipsoid (by minimizing maximum eigenvalue)
% subject to bound on worst case mode 2 cost: J < 1/epsilon
%epsilon = 0.1;

cvx_begin sdp
  cvx_solver mosek
  variable S(n, n) symmetric
  variables Y(m,n) bet(1)
  minimize(bet)
  subject to 
  CS_t = [S*C' zeros(size(C'))];

  for i=1:3
      for j = 1:4
          block1=[ ((A_term)*S+(B_term)*Y)', CS_t , Y';
            w(:,j)', zeros(1,m+n)];
          block2 = [S, Y';zeros(1,m+n)];
          block3 = blkdiag(S,bet);
          block4 = blkdiag(S,eye(size(Q)),inv(pp.R));
          
          [block3,block1;
           block1',block4] >= 0;
      end
  end
  for o=1:2
      block2=[(h(o))^2, F(o,:)*S+G(o,:)*Y; 
          (F(o,:)*S+G(o,:)*Y)', S];
      block2>=0; 
  end
cvx_end 

Qinv = S;
Q_N=inv(S);
K=Y*Q_N;

% cvx_begin sdp
%   cvx_solver mosek
%   variable gam(1)
%   minimize(-gam)
%   subject to 
%   for q = 1:2
%     F_ = (F(q,:)+G(q)*K);
% 
%     gam <= 1/(F_*Qinv*F_');
%   end
%   for i = 1:2
%     gam <= p.x_term^2/Qinv(i,i);
%   end
%   
%   gam*(Q+K'*R*K) >= double(bet)*Q_N;
% 
% cvx_end


gamma = inf;
for q = 1:2
  F_ = (F(q,:)+G(q)*K);
  gamma = min(gamma,h(q)^2/(F_*Qinv*F_'));
end
for i = 1:2
  gamma = min(gamma, p.x_term^2/Qinv(i,i));
end
Qisqrt = sqrtm(Qinv);
gamma_min = bet/max(eig(Qisqrt*(pp.Q+K'*pp.R*K)*Qisqrt));
if gamma < gamma_min
  error('Terminal constraint computation failed (not invariant: bet = %.4e beta_min = %.4e',bet,beta_min);
end
V = Q_N/gamma;
P = Q_N;

% V = Q_N/gam;
% P = Q_N;
