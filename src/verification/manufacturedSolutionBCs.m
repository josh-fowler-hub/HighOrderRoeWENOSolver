function [ULBC,URBC] = manufacturedSolutionBCs(U)
gamma = 1.4;
xs = zeros([1,3]);
ys = zeros([1,3]);
zs = zeros([1,3]);
L = 1;
W = manufacturedSolutionSetUp(L,xs,ys,zs);

ULBC = primToCons(gamma, W(1,:), W(2,:), W(5,:));
URBC = [U(:,end), U(:,end), U(:,end)];
end