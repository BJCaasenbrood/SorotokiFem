%------------------------------------ assemble global finite-element system 
function Fem = assembleGlobalFem(Fem,ForceBuild)
    
if nargin < 2
    ForceBuild = false; 
end

% % evaluate shape-functions at nodal locations
% if (Fem.Iteration == 1 && Fem.Increment == 1)
%     tab = struct; 
%     tab.Element = Fem.Element;
% 
%     switch(Fem.Mesh.Type)
%         case('C2PX'), tab = TabulateShapeFunctions(tab);
%         case('C2T3'), tab = TabulateShapeFunctions(tab);
%         case('C2Q4'), tab = TabulateShapeFunctions(tab);
%         %case('C3T3'), tab = TabulateShapeFunctions(tab);
%         case('C3H8'), tab = TabulateShapeFunctionsC3H8(tab);
%         case('C3T4'), tab = TabulateShapeFunctionsC3T4(tab);
%         otherwise,    tab = TabulateShapeFunctions(tab);
%     end
% 
%     %Fem.ShapeFnc = tab.ShapeFnc;
% end

dV  = 0*ones(Fem.Mesh.NElem,1);
% if Fem.VolumetricPressure
%     beta = Fem.OptFactor*Fem.LoadingFactor;
%     dV   = beta*Fem.PressureCell(:,2)*ones(Fem.NElem,1);
% elseif ~isempty(Fem.Contraction)
%     beta = Fem.OptFactor*Fem.LoadingFactor;
% 
%     if isa(Fem.Contraction,'function_handle')
%        dV = Fem.Contraction(Fem);
%     else
%        dV = beta*Fem.Contraction;
%     end
% end

beta    = Fem.options.loadingFactor;
[E,~,V] = materialFieldFem(Fem);

Fem.system.Potential = 0;
% Fem.PotentialG = 0;

index    = 0; 
subindex = 0;

% if isa(Fem.Material,'NeoHookean')
%     LocalType = 1;
% elseif strcmp(Fem.Material.Type,'Yeoh')
%     LocalType = 2;
% elseif strcmp(Fem.Material.Type,'Mooney')
%     LocalType = 3;
% end

if (~Fem.options.isAssembled && ~Fem.options.isNonlinear) ...
        || Fem.options.isNonlinear || ForceBuild
    
for el = 1:Fem.Mesh.NElem
   
    NDof = Fem.Mesh.ElemNDof(el);
    if Fem.Dim == 2
        eDof = reshape([2*Fem.Mesh.Element{el}-1;
               2*Fem.Mesh.Element{el}],NDof,1);
    else
        eDof = reshape([3*Fem.Mesh.Element{el}-2;
            3*Fem.Mesh.Element{el}-1;
            3*Fem.Mesh.Element{el}],NDof,1);
    end
 
    if isa(Fem.Material,'NeoHookean')
        nn = length(Fem.Mesh.Element{el});
        [Fe,Fb,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te,Ve,Vge,~,Re,Ue] = ...
            LocalsNHFast_mex(Fem.Mesh.Element{el},eDof,dV(el),full(E(el)),...
            Fem.Dim,Fem.Mesh.Node,...
            Fem.Mesh.ShapeFnc{nn}.N,...
            Fem.Mesh.ShapeFnc{nn}.dNdxi,...
            Fem.Mesh.ShapeFnc{nn}.W,...
            Fem.solver.sol.x, ...
            Fem.solver.sol.dx, ...
            Fem.Material.params.Rho, ...
            Fem.Material.params.Zeta,  ...
            Fem.system.Gravity,...
            Fem.Material.params.Mu,...
            Fem.Material.params.Lambda);
    elseif isa(Fem.Material,'Yeoh')
        nn = length(Fem.Element{el});
        [Fe,Fb,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te,Ve,Vge,~,Re,Ue] = ...
            LocalsYHFast(Fem.Element{el},eDof,dV(el),full(V(el)),...
            Fem.Dim,Fem.Node0,...
            Fem.ShapeFnc{nn}.N,...
            Fem.ShapeFnc{nn}.dNdxi,...
            Fem.ShapeFnc{nn}.W,...
            Fem.Utmp, Fem.dUtmp,...
            Fem.Material.Density, Fem.Material.Damping,...
            Fem.Gravity,...
            [Fem.Material.C1;Fem.Material.C2;Fem.Material.C3],...
            [Fem.Material.D1;Fem.Material.D2;Fem.Material.D3]);    
    elseif isa(Fem.Material,'Mooney')
        nn = length(Fem.Element{el});
        [Fe,Fb,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te,Ve,Vge,~,Re,Ue] = ...
            LocalsMNFast_mex(Fem.Element{el},eDof,dV(el),full(V(el)),...
            Fem.Dim,Fem.Node0,...
            Fem.ShapeFnc{nn}.N,...
            Fem.ShapeFnc{nn}.dNdxi,...
            Fem.ShapeFnc{nn}.W,...
            Fem.Utmp, Fem.dUtmp,...
            Fem.Material.Density, Fem.Material.MaterialDamping,...
            Fem.Gravity,...
            [Fem.Material.C10;Fem.Material.C01],...
            Fem.Material.K);       
     else
         [Fe,Fb,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te,Ve,Vge,~,Re,Ue] = ...
             Locals(Fem,Fem.Element{el},eDof,dV,V(el));
    end
    
%     if Fem.VolumetricPressure
%         dp = Fem.LoadingFactor*Fem.PressureCell(2)*invlerp(V(el),1,0);
%         [Fe2,~,~,~,~,Ke2,Kte2] = ...
%             Locals(Fem,Fem.Element{el},eDof,dp,V(el));
%         Fb  = Fb + Fe2/E(el);
%         Ke  = Ke + Ke2/E(el);
%         Kte = Kte + Kte2/E(el);
%     end
     
    ind1 = index+1:index+NDof^2;               % matrix indexing
    ind2 = index+1:index+NDof;                 % vector indexing
    ind3 = subindex+1:subindex+NDof/Fem.Dim;   % dimension indexing

    I = repmat(eDof,1,NDof); J = I';
    Fem.triplets.e(ind1)  = el;
    Fem.triplets.i(ind1)  = I(:);
    Fem.triplets.j(ind1)  = J(:);
    Fem.triplets.m(ind1)  = Me(:);
    Fem.triplets.c(ind1)  = Ce(:);
    Fem.triplets.k(ind1)  = Ke(:);
    Fem.triplets.t(ind1)  = Kte(:);
    Fem.triplets.fi(ind2) = Fe(:);
    Fem.triplets.fb(ind2) = Fb(:);
    Fem.triplets.ft(ind2) = Te(:);    
    
    Fem.triplets.s(ind3,1) = E(el)*Svme(:);
    Fem.triplets.s(ind3,2) = E(el)*SS(:,1);
    Fem.triplets.s(ind3,3) = E(el)*SS(:,2);
    Fem.triplets.s(ind3,4) = E(el)*SS(:,4);
    Fem.triplets.p(ind3,1) = E(el)*EE(:,1);
    Fem.triplets.p(ind3,2) = E(el)*EE(:,2);
    Fem.triplets.p(ind3,3) = E(el)*EE(:,4);
    Fem.triplets.l(ind3)   = Fem.Mesh.Element{el}(:);
    
    Fem.system.Potential  = Fem.system.Potential + beta*Ve;
    %Fem.PotentialG = Fem.PotentialG + beta*Vge;
    Fem.triplets.v(ind3) = Qe(:);
    
    if Fem.solver.Iteration == 1
        Fem.triplets.ElemRot{el} = Re;
        Fem.triplets.ElemStr{el} = Ue;
    end
    
    index    = index + NDof^2;
    subindex = subindex + NDof/Fem.Dim;
end

end

qa = Fem.system.Ia;

%M  = sparse(Fem.i,Fem.j,E(Fem.e).*Fem.m);  % init mass matrix
Fem.system.Mass = sparse(Fem.triplets.i,Fem.triplets.j,Fem.triplets.m);
Fem.system.Mass = Fem.system.Mass(qa,qa);
% init mass matrix
Fem.system.Damping = sparse(Fem.triplets.i,Fem.triplets.j,Fem.triplets.c);  
Fem.system.Damping = Fem.system.Damping(qa,qa);


Fem.system.Stiffness = sparse(Fem.triplets.i,Fem.triplets.j,           ...
    E(Fem.triplets.e).*Fem.triplets.k);  % init global stiffness
Fem.system.Stiffness = Fem.system.Stiffness(qa,qa);

Fem.system.Tangent   = sparse(Fem.triplets.i,Fem.triplets.j,           ...
    E(Fem.triplets.e).*Fem.triplets.t);  % init global tangent stiffness
Fem.system.Tangent = Fem.system.Tangent(qa,qa);

Fem.system.fElastic  = sparse(Fem.triplets.i,1,E(Fem.triplets.e)       ...
    .*Fem.triplets.fi);     % init elastic force vector

Fem.system.fElastic = Fem.system.fElastic(qa);

Fem.system.fBody     = sparse(Fem.triplets.i,1,E(Fem.triplets.e)       ...
    .*Fem.triplets.fb);     % init body force vector
Fem.system.fBody = beta * Fem.system.fBody(qa);


Fem.options.isAssembled = true;

end