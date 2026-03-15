function [Ws, inits] = manufactured_solution_set_up(L, x, y, z)
W = manufactured_solution(2);
rs = W(1,:);
us = W(2,:);
vs = W(3,:);
ws = W(4,:);
ps = W(5,:);
r0 = ones(size(x)).*rs(1);
u0 = ones(size(x)).*us(1);
v0 = ones(size(x)).*vs(1);
w0 = ones(size(x)).*ws(1);
p0 = ones(size(x)).*ps(1);


% r = rs(1) + rs(2) .* sin((rs(5) .* pi .* x)./L) + rs(3) .* cos((rs(6) .* pi .* y)./L) + rs(4) .* sin((rs(7) .* pi .* z)./L);
% 
% u = us(1) + us(2) .* sin((us(5) .* pi .* x)./L) + us(3) .* cos((us(6) .* pi .* y)./L) + us(4) .* cos((us(7) .* pi .* z)./L);
% 
% v = vs(1) + vs(2) .* cos((vs(5) .* pi .* x)./L) + vs(3) .* sin((vs(6) .* pi .* y)./L) + vs(4) .* sin((vs(7) .* pi .* z)./L);
% 
% w = ws(1) + ws(2) .* sin((ws(5) .* pi .* x)./L) + ws(3) .* sin((ws(6) .* pi .* y)./L) + ws(4) .* cos((ws(7) .* pi .* z)./L);
% 
% p = ps(1) + ps(2) .* cos((ps(5) .* pi .* x)./L) + ps(3) .* sin((ps(6) .* pi .* y)./L) + ps(4) .* cos((ps(7) .* pi .* z)./L);


r = rs(1) + rs(2) .* sin((rs(5) .* pi .* x)./L);

u = us(1) + us(2) .* sin((us(5) .* pi .* x)./L);

v = vs(1) + vs(2) .* cos((vs(5) .* pi .* x)./L);

w = ws(1) + ws(2) .* sin((ws(5) .* pi .* x)./L);

p = ps(1) + ps(2) .* cos((ps(5) .* pi .* x)./L);

Ws = [r; u; v; w; p];
inits = [r0; u0; v0; w0; p0];

end