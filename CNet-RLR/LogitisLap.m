function theta = LogitisLap (X,y,L,lambda1,lambda2)
% global p 
%% To get the theta which minimize the cost function -sum_i{y_i*log(g(theta'*x_i))+(1-y_i)*log(1-g(theta'*x_i))}
%% +lambda1 |theta|_1+lambda2 theta' A theta
%% The input y is a n by 1 matrix
%% Input X is G X n gene expression matrix
%% L is graphical Laplacian matrix, lambda1 and lambda2 are parameters in the cost function that are known
%% The output theta is the optimal solution of the covex optimization problem
%% The dimension of theta: G by 1
[dimn,~]=size(y); % get the number of n and set it to dimn
[dimG,m]=size(X); % get the number of G and set it to dimG

%% 
% options = optimoptions('fmincon','Algorithm','interior-point');
% 
% initial_theta = zeros(p,1);
% A = [];
% b = [];
% Aeq = []; 
% beq = [];
% vlb = -1*ones(1,p); 
% vub =  1*ones(1,p);
% 
% [theta, cost]  = ...
%     fmincon(@(theta)(myfunction(theta, X, y, lambda1, lambda2)), initial_theta, A, b, Aeq, beq, vlb, vub, [], options);
% 


%%  To use cvx to solve the convex problem
cvx_begin
    variables theta(dimG) ;

    minimize (  sum((1-y)'.*(theta'*X))+sum( log_sum_exp( [zeros(1,m); -theta'*X]) ) +lambda1*norm(theta,1)+lambda2*theta'*L*theta  );
%     minimize (  sum((1-y)'.*(theta'*X))+sum( log_sum_exp( -theta'*X ) ) +lambda1*norm(theta,1)+lambda2*theta'*L*theta  );
%       minimize (  sum((1-y)'.*(theta'*X))+sum( log_sum_exp( -theta'*X ) ) +lambda1*norm(theta,1)+lambda2*theta'*L*theta  );
cvx_end
    
    
