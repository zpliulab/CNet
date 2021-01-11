function J = costFunction12(theta, u)
global  L1 X_train  y_train lambda_opt alpha_opt 


% m = length(y); % number of training examples

z = X_train * theta;
hx = 1 ./ (1 + exp(-z));

R1 = sum(u);
% R1 = sum(abs(theta));
R2 = theta'*L1*theta/2;

J = sum([- y_train' * log(hx) - (1 -  y_train)' * log(1 - hx)]) + lambda_opt*alpha_opt*R1 + lambda_opt*(1-alpha_opt)*R2;

return