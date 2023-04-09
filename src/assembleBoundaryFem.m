function Fem = assembleBoundaryFem(Fem)

Ft  = sparse(Fem.Mesh.NNode*Fem.Dim,1);  % init contraction force
Fnc = sparse(Fem.Mesh.NNode*Fem.Dim,1);  % init normal contact force
Ftc = sparse(Fem.Mesh.NNode*Fem.Dim,1);  % init tangent contact force
F   = sparse(Fem.Mesh.NNode*Fem.Dim,1);  % init global force vector
L   = sparse(Fem.Dim*Fem.Mesh.NNode,1);  % init output vector

if Fem.options.isNonlinear && ~isempty(Fem.options.loadingFactor)
    %beta = Fem.OptFactor*Fem.options.loadingFactor;
    beta = Fem.options.loadingFactor;
else
    beta = 1;
end

% adding basic loads
if isfield(Fem.system,'Load') && ~Fem.options.isPrescribed
    F = beta * addBasicLoadsFem(Fem,F);
end

%%
qa = Fem.system.Ia;
Fem.system.fContact   = Fnc(qa);
Fem.system.fTangent   = Ftc(qa);
%Fem.system.fInternal  = Fe; 
Fem.system.fInput     = F(qa) + Ft(qa) + Fnc(qa) + Ftc(qa);
% Fem.Stiffness        = K + spMat;
% Fem.TangentStiffness = Ktr + spMat;
%Fem.Residual         = Fem.fInternal - Fem.fExternal - Fem.fInput;
Fem.system.Output    = L;

end

