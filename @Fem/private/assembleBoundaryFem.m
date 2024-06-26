function Fem = assembleBoundaryFem(Fem)

qa  = Fem.system.Ia;
Ft  = sparse(Fem.Mesh.NNode*Fem.Dim,1);  % init contraction force
Fnc = sparse(Fem.Mesh.NNode*Fem.Dim,1);  % init normal contact force
Ftc = sparse(Fem.Mesh.NNode*Fem.Dim,1);  % init tangent contact force
Fc  = sparse(Fem.Mesh.NNode*Fem.Dim,1);  % init global force vector
F   = sparse(Fem.Mesh.NNode*Fem.Dim,1);  % init global force vector
L   = sparse(Fem.Dim*Fem.Mesh.NNode,1);  % init output vector

beta = 1.0; % loading increment parameter
if Fem.options.isNonlinear && ~isempty(Fem.options.loadingFactor)
    beta = Fem.options.loadingFactor;
end

% adding basic loads
if isfield(Fem.system,'Load') && ~Fem.options.isPrescribed
    F = beta * assembleLoadsFem(Fem);
end

% adding basic loads
if isfield(Fem.system,'Tendon') && ~Fem.options.isPrescribed
    F = beta * assembleTendonForces(Fem);
end

% adding basic loads
if isfield(Fem.system,'Pressure') && ~Fem.options.isPrescribed
    F = beta * assemblePressureForces(Fem);
end

% adding contact loads
if isfield(Fem.system,'Contact') && ~Fem.options.isPrescribed
    [Fnc, Ftc, Knc, Ktc] = assembleContactForces(Fem);
    Fem.system.Tangent = Fem.system.Tangent ... 
        + Knc(qa,qa) + Ktc(qa,qa);

    if Fem.options.isForceContactDamping        
        Fem.system.Damping = Fem.system.Damping ... 
            + Knc(qa,qa) * Fem.solver.TimeStep;     
    end 

end

% adding displacement loads
if isfield(Fem.system,'Displace') && ~Fem.options.isPrescribed
    [ge, ce] = assembleDisplaceEq(Fem);
    cDofs   = find(ce == 1);
    cMatrix = diag(ce);
    AllDofs = 1:Fem.Dim*Fem.Mesh.NNode;

    [Ice,~]  = ismember(AllDofs(:),cDofs(:));

    Fem.system.cMatrix   = cMatrix(cDofs,qa);
    Fem.system.cResidual = ge(cDofs);
    Fem.system.Ic        = Ice;
end

% adding controller loads
if isfield(Fem.system,'Controller') 
    u = Fem.system.Controller(Fem);
    G = Fem.system.InputMap(Fem);
    Fc = G * u;

    Fem.system.fControl = u;
elseif isfield(Fem.system,'fControl') && ~isfield(Fem.system,'Controller') 
    G  = Fem.system.InputMap(Fem);
    Fc = G * Fem.system.fControl;
end

% adding output
if isfield(Fem.system,'Output') && ~Fem.options.isPrescribed
    L = assembleOutputFem(Fem);
    kdof = find(L ~= 0);
    fdof = find(F ~= 0);
    ok = find(Fem.triplets.i == kdof);
    fk = find(Fem.triplets.i == fdof);
    Fem.triplets.k(ok) = Fem.topology.OutputStiffness;
    Fem.triplets.k(fk) = Fem.topology.InputStiffness;
end

if  isfield(Fem.system,'Spring') && ~Fem.options.isPrescribed
    [Fs, Ks] = assembleSpringForcesFem(Fem);
    F = Fs + F;
    Fem.system.Tangent = Fem.system.Tangent + Ks(qa,qa);
end

Fem.system.fContact = Fnc(qa);
Fem.system.fTangent = Ftc(qa);
Fem.system.fInput   = F(qa) + Ft(qa) + Fnc(qa) + Ftc(qa) + Fc(qa);
Fem.system.fOutput  = L;
end

