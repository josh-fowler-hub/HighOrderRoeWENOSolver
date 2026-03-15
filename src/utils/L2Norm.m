function [L2] = L2Norm(U, Ustar, NN)
for i = 1:length(U(1,:))
    L21 = sqrt((sum((U(1,:) - Ustar(1,:)).^2))/NN);
    L22 = sqrt((sum((U(2,:) - Ustar(2,:)).^2))/NN);
    L23 = sqrt((sum((U(3,:) - Ustar(3,:)).^2))/NN);
end
L2 = [L21; L22; L23];
end