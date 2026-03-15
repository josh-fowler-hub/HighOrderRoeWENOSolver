function dt = calculateTimeStep(cfl, dx, time, tfinal, U, gamma)
% CFL-based time step calculation
smax = findWaveSpeed(U, gamma);

dt = (cfl*dx)/smax;
if isnan(dt)
    dt = 1e-6;
end

if time + dt > tfinal
    dt = tfinal - time;
end
end
