function Fem = assembleBoundaryFem(Fem)

qa  = Fem.system.Ia;
Ft  = sparse(Fem.Mesh.NNode*Fem.Dim,1);  % init contraction force
Fnc = sparse(Fem.Mesh.NNode*Fem.Dim,1);  % init normal contact force
Ftc = sparse(Fem.Mesh.NNode*Fem.Dim,1);  % init tangent contact force
F   = sparse(Fem.Mesh.NNode*Fem.Dim,1);  % init global force vector
L   = sparse(Fem.Dim*Fem.Mesh.NNode,1);  % init output vector

if Fem.options.isNonlinear && ~isempty(Fem.options.loadingFactor)
    beta = Fem.options.loadingFactor;
else
    beta = 1;
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
end

% adding displacement loads
if isfield(Fem.system,'Displace') && ~Fem.options.isPrescribed
    [ge, ce] = assembleDisplaceEq(Fem);
    cDofs   = find(ce == 1);
    cMatrix = diag(ce);
    AllDofs    = 1:Fem.Dim*Fem.Mesh.NNode;

    [Ice,~]  = ismember(AllDofs(:),cDofs(:));

    Fem.system.cMatrix   = cMatrix(cDofs,qa);
    Fem.system.cResidual = ge(cDofs);
    Fem.system.Ic        = Ice;
end

%
Fem.system.fContact = Fnc(qa);
Fem.system.fTangent = Ftc(qa);
Fem.system.fInput   = F(qa) + Ft(qa) + Fnc(qa) + Ftc(qa);
% Fem.Stiffness        = K + spMat;
% Fem.TangentStiffness = Ktr + spMat;
%Fem.fResidual     = Fem.fInternal - Fem.fExternal - Fem.fInput;
Fem.system.Output = L;
end

