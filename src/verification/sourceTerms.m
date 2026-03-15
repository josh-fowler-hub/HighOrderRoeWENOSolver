function fx = sourceTerms(x)
W = manufacturedSolution(2);
rs = W(1,:);
us = W(2,:);
ps = W(5,:);
L = 1;

r0 = rs(1);
rx = rs(2);
arx = rs(5);
u0 = us(1);
ux = us(2);
aux = us(5);
p0 = s(1);
px = ps(2);
apx = ps(5);
gamma = 1.4;



fx(1,:) = (aux.*ux.*pi.*cos((pi.*aux.*x)./L).*(r0 + rx.*sin((pi.*arx.*x)./L)))./L + (arx.*rx.*pi.*cos((pi.*arx.*x)./L).*(u0 + ux.*sin((pi.*aux.*x)./L)))./L;
fx(2,:) = (arx.*rx.*pi.*cos((pi.*arx.*x)./L).*(u0 + ux.*sin((pi.*aux.*x)./L)).^2)./L - (apx.*px.*pi.*sin((pi.*apx.*x)./L))./L + (2.*aux.*ux.*pi.*cos((pi.*aux.*x)./L).*(r0 + rx.*sin((pi.*arx.*x)./L)).*(u0 + ux.*sin((pi.*aux.*x)./L)))./L;
fx(3,:) = (aux.*ux.*pi.*cos((pi.*aux.*x)./L).*(p0 + px.*cos((pi.*apx.*x)./L) + (r0 + rx.*sin((pi.*arx.*x)./L)).*((u0 + ux.*sin((pi.*aux.*x)./L)).^2./2 + (p0 + px.*cos((pi.*apx.*x)./L))./((r0 + rx.*sin((pi.*arx.*x)./L)).*(gamma - 1)))))./L - (u0 + ux.*sin((pi.*aux.*x)./L)).*((r0 + rx.*sin((pi.*arx.*x)./L)).*((apx.*px.*pi.*sin((pi.*apx.*x)./L))./(L.*(r0 + rx.*sin((pi.*arx.*x)./L)).*(gamma - 1)) - (aux.*ux.*pi.*cos((pi.*aux.*x)./L).*(u0 + ux.*sin((pi.*aux.*x)./L)))./L + (arx.*rx.*pi.*cos((pi.*arx.*x)./L).*(p0 + px.*cos((pi.*apx.*x)./L)))./(L.*(r0 + rx.*sin((arx.*x.*pi)./L)).^2.*(gamma - 1))) + (apx.*px.*pi.*sin((pi.*apx.*x)./L))./L - (arx.*rx.*pi.*cos((pi.*arx.*x)./L).*((u0 + ux.*sin((pi.*aux.*x)./L)).^2./2 + (p0 + px.*cos((pi.*apx.*x)./L))./((r0 + rx.*sin((pi.*arx.*x)./L)).*(gamma - 1))))./L);

end