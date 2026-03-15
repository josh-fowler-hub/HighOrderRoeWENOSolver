syms x r0 rx arx u0 ux aux p0 px apx L gamma t;

r = r0 + rx * sin((arx * pi * x)/L);

u = u0 + ux * sin((aux * pi * x)/L);

p = p0 + px * cos((apx * pi * x)/L);

es = p/(r.*(gamma - 1));
U = [r; r*u; r*(0.5*u^2 + es)];
F = [r*u; r*u^2 + p; u*(r*(0.5*u^2 + p/((gamma-1)*r)) + p)];

fx = diff(U,t) + diff(F,x)