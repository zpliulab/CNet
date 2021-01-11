function [ AUC, pp, y] = Predict(train_data, train_y, Est_theta)

X = train_data; 
y0 = train_y; 
theta = Est_theta;
[n,p] = size(X);
% if(!is.matrix(w)){w = matrix(c(w),ncol = 1)}
% X = [X, ones(n,1)];
pp = 1./(1+exp(-X*theta));
y = ones(n,1);
for i = 1:n
    if pp(i) < 0.5
        y(i) = 0;
    end
end
% acc = length(which(y==y0))/n; 
% AUC = accuracy.rate
% AUC = auc(y0, c(pp))[1];
AUC = plot_roc(pp,y0);

return