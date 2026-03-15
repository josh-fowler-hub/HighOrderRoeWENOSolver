function [coeffs] = SSPRKcoeffs(order)
stages = order + 1;
zero = 1;
alpha = zeros([order,zero+stages-1]);
alpha(2,1) = 0;
alpha(2,2) = 1;

for m = 3:order
    for k = 1:m-2
        alpha(m,k+1) = (2/k)*alpha(m-1,k);
    end
    alpha(m,m) = (2/m)*alpha(m-1,m-1);
    alpha(m,1) = 1 - sum(alpha(m,2:m));
end
coeffs = alpha(order,:);
end