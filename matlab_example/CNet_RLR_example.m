clc; clear; close all; format short;  % format long;
global L1  X_train  y_train  lambda_opt alpha_opt 
%% 邻接矩阵
adj = importdata('C:\Users\LiLingyu\Desktop\NSLR-master\adj_example.csv');
A_adj = sparse(adj.data);
%% Set folds in c-v
nfold = 5;
%% provide the value of alpha
% alpha = 1;     % --- lasso
alpha = 0.3; %    1e-5;   % --- lG
% L= eye(p,p);
% alpha = 0.5;   % -- elastic net
%% 割点
vector_hat = importdata('C:\Users\LiLingyu\Desktop\NSLR-master\vector_hat_example.txt'); 
delta_hat = vector_hat.data;
%% Set options
options = optimoptions('fmincon','Algorithm','interior-point');
options.MaxFunEvals = 1e5;
%% 循环 30 次
% for i = 1
i = 1; 
%% 数据输入
txt = importdata(['C:\Users\LiLingyu\Desktop\NSLR-master\Data_train\',num2str(i),'.txt']);  
train_data = txt.data;
[m,q] = size(train_data);
X_train = train_data(1:m,2:q);     % 以索引的前1000个数据点作为测试样本Xtest
y_train = train_data(1:m,1);  
p = q-1;
%% laplace 矩阵
[L, L1] = Laplacian_Matrix(p,A_adj);
% L = eye(p,p);
%% get the least upper bound, lambda_max
% lammax = getLambMax(X_train, y_train, alpha); 
% e = (log(lammax)-log(1))/19;
% lambda = exp(log(1):e:log(lammax)); 
%% get the optimal lambda and alpha through cross validation 
% [lambda_opt, alpha_opt, r] = cv( X_train', y_train, L1, alpha, lambda, nfold );
lambda_opt = 60;
%% 最优参数
% alpha_opt = 0.5;
% lambda_opt = 4.874056378116327e+02;   % L1
%% 约束条件
theta_0 = zeros(p,1); % 初值为0
u_0 = (1e-4)*ones(p,1);  % 初值为1

a11 = eye(p,p);
a22 = -1*eye(p,p);
A1 = [a11;a22];
A2 = [a22;a22];
A = [A1,A2];
b = zeros(2*p,1);

% A = [];
% b = [];
Aeq = []; 
beq = [];
% vlb = [];        
% vub = [];  

vlb1 = zeros(p,1);        
vub1 = zeros(p,1);  
delta1 = abs(delta_hat);
for j = 1:p
    if delta1(j) == 1
        vlb1(j) = 1e-5;
        vub1(j) = 1;
    else 
        vlb1(j) = -inf;
        vub1(j) =  inf;
    end
end
vlb = [vlb1; vlb1];        
vub = [vub1; vub1];  
%% 内点法求解
[X_sol, cost, exitflag, output, mu, grad, hessian] = fmincon(@(x)(costFunction12(x(1:p),x(p+1:2*p))), [theta_0;u_0], A, b, Aeq, beq, vlb, vub, [], options);

theta = X_sol(1:p);
u = X_sol(p+1:2*p);
%% 输出 Cost 和 theta
fprintf('Cost at theta found by fmincon: %e', cost);
fprintf('\n')
%% 预测
txt1 = importdata(['C:\Users\LiLingyu\Desktop\NSLR-master\Data_test\',num2str(i),'.txt']);   
test_data = txt1.data;
[k,l] = size(test_data);
X_test = test_data(1:k,2:l);   
y_test = test_data(1:k,1);  

[AUC_test, y_pre, y_true] = Predict(X_test, y_test, theta); 
% [AUC_test, y_pre, y_true] = Predict(X_train, y_train, theta);  
P_test = [y_pre, y_true];
para = [lambda_opt, AUC_test, cost]
%% 保存
% csvwrite(['C:\Users\LiLingyu\Desktop\NSLR-master\theta',num2str(i),'.csv'],theta);
% csvwrite(['C:\Users\LiLingyu\Desktop\NSLR-master\P_test',num2str(i),'.csv'],P_test);
% csvwrite(['C:\Users\LiLingyu\Desktop\NSLR-master\para',num2str(i),'.csv'],para);
% end