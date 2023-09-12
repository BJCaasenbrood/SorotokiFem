function dz = flow(Fem, t, z, varargin)

    if ~isfield(Fem.system,'Ia')
        qa = 1:Fem.Dim * Fem.Mesh.NNode;
    else
        qa = Fem.system.Ia;
    end
    
    nq = numel(z) / 2;
    
    Fem.solver.Time      = t;
    Fem.solver.sol.x(qa)  = z(1:nq);
    Fem.solver.sol.dx(qa) = z(nq+1:2*nq);

    Fem = Fem.compute();

    dz  = zeros(2 * nq, 1);
    dz(1:nq)      = z(nq+1:2*nq);
    dz(nq+1:2*nq) = - Fem.system.Mass \ Fem.system.fResidual;

    cla;
    showVonMisesFem(Fem);
    drawnow;
    t

end