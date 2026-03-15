function U = primToCons(gamma, rhos, us, Ps)
es = Ps./(rhos.*(gamma - 1));
U = [rhos; rhos.*us; rhos.*(0.5.*us.^2 + es)];
end
