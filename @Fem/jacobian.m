function J = jacobian(Fem, t, z, varargin)

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

    J  = zeros(2 * nq, 2 * nq);
    J(1:nq,nq+1:2*nq) = eye(nq);
    J(nq+1:2*nq,1:nq) = -Fem.system.Mass \ Fem.system.Tangent;
    J(nq+1:2*nq,nq+1:2*nq) = -Fem.system.Mass \ (Fem.system.Damping);
end