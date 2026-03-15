function F = Fvec(r,u,P,g)
F(1,:) = r.*u;
F(2,:) = r.*u.^2+P;
F(3,:) = u.*(r.*(0.5.*u.^2+P./((g-1).*r))+P);
end