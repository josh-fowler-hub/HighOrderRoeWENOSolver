function F_vec = dUdt(U, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time)
% dUdt: Computes the flux vector to be added to the conservative variable vector.
% Inputs:
%   U: conservative variables (3xN)
%   flux: flux limiter identifier
%   lam: scaling factor for MUSCL
%   gamma: ratio of specific heats
%   dx: spatial step
%   kap: MUSCL parameter
%   rec: reconstruction method
%   CD_Term_Order: central differencing order
%   riemann_solver: type of Riemann solver
%   xs: x-grid vector
%   time: current time (some solvers or BCs may use it)
% Output:
%   F_vec: flux vector (3xN)

switch rec
    case 'ROE'
        [URi,ULi,URim1,ULim1] = musclRec(U,18,lam,kap);
        BCs = [0 0 0; 0 0 0; 0 0 0];
    case 'MUSCL'
        [URi,ULi,URim1,ULim1] = musclRec(U,flux,lam,kap);
        BCs = [0 0 0; 0 0 0; 0 0 0];
    case 'WENO3'
        [URi,ULi] = weno3Rec(U(:,2:end));
        [URim1,ULim1] = weno3Rec(U(:,1:end-1));
        BCs = [0 0 0; 0 0 0; 0 0 0];
    case 'WENO5'
        [URi,ULi] = weno5Rec(U(:,2:end));
        [URim1,ULim1] = weno5Rec(U(:,1:end-1));
        BCs = [0 0 0; 0 0 0; 0 0 0];
    otherwise
        error('Unknown reconstruction method: %s', rec);
end

switch riemann_solver
    case 'Roe'
        [Fip12, Fim12] = roeSolver(U,URi,ULi,URim1,ULim1,gamma,CD_Term_Order);
    case 'Roe-Pike'
        [Fip12, Fim12] = roePikeSolver(U,URi,ULi,URim1,ULim1,gamma,CD_Term_Order);
    case 'HLL'
        [Fip12, Fim12] = hllSolver(U,URi,ULi,URim1,ULim1,xs,gamma,time);
    case 'HLLC'
        [Fip12, Fim12] = hllcSolver(U,URi,ULi,URim1,ULim1,xs,gamma,time);
    case 'Osher'
        [Fip12, Fim12] = osherSolver(U,URi,ULi,URim1,ULim1,xs,gamma,time);
    otherwise
        error('Unknown Riemann solver: %s', riemann_solver);
end

F_vec = -(1/dx)*(Fip12 - Fim12);
F_vec = [BCs, F_vec, BCs];
end
