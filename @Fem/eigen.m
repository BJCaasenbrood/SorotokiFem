function Fem = eigen(Fem,varargin)
    if isempty(varargin)
        x = Fem.solver.sol.x;
    else
        x = varargin{1};
    end

    % precompute
    Fem = Fem.compute(x);

    % get 
    ia = Fem.system.Ia; 
    M = Fem.system.Mass;
    K = Fem.system.Tangent;

    [V, D] = eigs(M,K, 10);
    Z = zeros(Fem.Mesh.NNode * Fem.Dim , 10);
    % Fem.solver.sol.yout       = Z;
    % Fem.solver.sol.yout(ia,:) = V;
    Fem.solver.sol.pod.V = Z;
    Fem.solver.sol.pod.V(ia,:) = V;
    Fem.solver.sol.pod.D = diag(D);
end