function U = SSPRKN(order,U,dUdt,dt)
alpha = SSPRKcoeffs(order);
Uvec = [U];
for i = 1:order
    Ui = Uvec(i,:) + dt/2 * dUdt(Uvec(i,:));
    Uvec = [Uvec;Ui];
end
Um = zeros(size(U));
for i = 1:length(alpha)-1
    Um = Um + alpha(i)*Uvec(i,:);
end

Um = Um + alpha(end-1)*U;

end