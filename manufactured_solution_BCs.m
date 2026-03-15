function [ULBC,URBC] = manufactured_solution_BCs(U)
gamma = 1.4;
xs = zeros([1,3]);
ys = zeros([1,3]);
zs = zeros([1,3]);
L = 1;
W = manufactured_solution_set_up(L,xs,ys,zs);

ULBC = prim_to_cons(gamma, W(1,:), W(2,:), W(5,:));
URBC = [U(:,end), U(:,end), U(:,end)];
end