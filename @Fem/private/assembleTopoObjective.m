function [f,dfdE,dfdV] = assembleTopoObjective(Fem)

X      = Fem.solver.sol.x;
fDof   = Fem.system.FreeDofs;
qa     = Fem.system.Ia;
E = materialFieldFem(Fem);
% [E,dEdy,~,dVdy] = materialFieldFem(Fem);

if strcmpi(Fem.topology.Type,'compliance') &&  ~Fem.options.isNonlinear

    f = (Fem.system.fInput + Fem.system.fBody).'*X(qa);
    temp = cumsum(-X(Fem.triplets.i).*Fem.triplets.k.*X(Fem.triplets.j));
    temp = temp(cumsum(Fem.Mesh.ElemNDof.^2));
    dfdE = [temp(1);temp(2:end)-temp(1:end-1)];

elseif strcmpi(Fem.topology.Type,'compliance') &&  Fem.options.isNonlinear

    K = sparse(Fem.triplets.i,Fem.triplets.j, ...
        E(Fem.triplets.e).*Fem.triplets.t);

    f = (Fem.system.fInput + Fem.system.fBody).'*X(qa);
    lam = 0*X;
    lam(qa) = K(qa,qa)\Fem.system.fInput;
    temp = cumsum(-X(Fem.triplets.i).*Fem.triplets.k.* ...
        lam(Fem.triplets.j));
    temp = temp(cumsum(Fem.Mesh.ElemNDof.^2));
    dfdE = [temp(1);temp(2:end)-temp(1:end-1)];

elseif strcmpi(Fem.topology.Type,'compliant') && ~Fem.options.isNonlinear

    L = Fem.system.fOutput;
    f = -L.'*X / (numel(L(abs(L)>0)));
    K = sparse(Fem.triplets.i,Fem.triplets.j, E(Fem.triplets.e) ...
        .*Fem.triplets.k);

    lam = 0 * X;
    lam(fDof) = K(fDof,fDof) \ L(fDof);
    temp = cumsum(X(Fem.triplets.i,1).*Fem.triplets.k.*lam(Fem.triplets.j,1));
    temp = temp(cumsum(Fem.Mesh.ElemNDof.^2));
    dfdE = [temp(1);temp(2:end)-temp(1:end-1)];

elseif strcmpi(Fem.topology.Type,'compliant') && Fem.options.isNonlinear

    L = Fem.OutputVector;
    f = -L.'* X / (numel(L(abs(L)>0)));
    K = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.t);
    lam = 0 * X;
    lam(fDof) = K(fDof,fDof)\Fem.OutputVector(fDof);
    temp = cumsum(u(Fem.i).*Fem.k.*lam(Fem.j));
    temp = temp(cumsum(Fem.ElemNDof.^2));
    dfdE = [temp(1);temp(2:end)-temp(1:end-1)];

end

% if Fem.VolumetricPressure
%     temp = cumsum(Fem.dFdE(Fem.e).*lam(Fem.j));
%     temp = temp(cumsum(Fem.ElemNDof.^2));
%     dfdE = dfdE + [temp(1);temp(2:end)-temp(1:end-1)];
% end

dfdV = zeros(size(Fem.Mesh.NElem));

end