function [f,dfdE,dfdV] = assembleObjectiveFem(Fem)

X      = Fem.solver.sol.x;
fDof   = Fem.system.FreeDofs;
qa     = Fem.system.Ia;
[E,dEdy,~,dVdy] = materialFieldFem(Fem);
type  = Fem.topology.Type;

if strcmpi(type,'compliance') &&  ~Fem.options.isNonlinear

    f = (Fem.system.fInput + Fem.system.fBody).'*X(qa);
    temp = cumsum(-X(Fem.triplets.i).*Fem.triplets.k.*X(Fem.triplets.j));
    temp = temp(cumsum(Fem.Mesh.ElemNDof.^2));
    dfdE = [temp(1);temp(2:end)-temp(1:end-1)];

elseif strcmpi(type,'compliance') &&  Fem.options.isNonlinear

    K = sparse(Fem.triplets.i,Fem.triplets.j, ...
        E(Fem.triplets.e).*Fem.triplets.t);
    f = (Fem.system.fInput + Fem.system.fBody).'*X(qa);
    lam = 0*X;
    lam(qa) = K(qa,qa)\Fem.system.fInput;
    temp = cumsum(-X(Fem.triplets.i).*Fem.triplets.k.* ...
        lam(Fem.triplets.j));
    temp = temp(cumsum(Fem.Mesh.ElemNDof.^2));
    dfdE = [temp(1);temp(2:end)-temp(1:end-1)];

elseif strcmpi(type,'compliant') && ~Fem.options.isNonlinear

    E = MaterialField(Fem);
    L = Fem.OutputVector;
    f = -L.'*u(:,1)/(numel(L(abs(L)>0)));
    K = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.t);
    lam = 0*u(:,1);
    lam(fDof) = K(fDof,fDof)\Fem.OutputVector(fDof);
    temp = cumsum(u(Fem.i,1).*Fem.k.*lam(Fem.j,1));
    temp = temp(cumsum(Fem.ElemNDof.^2));
    dfdE = [temp(1);temp(2:end)-temp(1:end-1)];

elseif strcmpi(type,'compliant') && Fem.options.isNonlinear

    E = MaterialField(Fem);
    L = Fem.OutputVector;
    f = -L.'*u(:,1)/(numel(L(abs(L)>0)));
    K = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.t);
    lam = 0*u(:,1);
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