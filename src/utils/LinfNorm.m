function [Linf] = LinfNorm(U, Ustar)

for i = 1:length(U(1,:))
    Linf1 = max(abs(U(1,:) - Ustar(1,:)));
    Linf2 = max(abs(U(2,:) - Ustar(2,:)));
    Linf3 = max(abs(U(3,:) - Ustar(3,:)));
end
Linf = [Linf1; Linf2; Linf3];
end