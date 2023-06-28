%------------------------------------ assemble global finite-element system 
function Fem = assembleGlobalFem(Fem,ForceBuild)
    
if nargin < 2
    ForceBuild = false; 
end

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

if isfield(Fem.system,'Dilation')
    for kk = 1:size(Fem.system.Dilation,1)
        elDof = Fem.system.Dilation{kk,1};
        eldV  = Fem.system.Dilation{kk,2}(Fem.solver.Time);

        dV(elDof) = eldV;
    end
end

beta    = Fem.options.loadingFactor;
[E,~,V] = materialFieldFem(Fem);

Fem.system.Potential = 0;

index    = 0; 
subindex = 0;

if (~Fem.options.isAssembled && ~Fem.options.isNonlinear) || Fem.options.isNonlinear

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

    EMat = Fem.materials.MatElem(el);
    Material = Fem.materials.Material{EMat};
 
    if isa(Material,'NeoHookean')
        nn = length(Fem.Mesh.Element{el});
        [Fe,Fb,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te,Ve,Vge,~,Re,Ue] = ...
            LocalsNHFast_mex(Fem.Mesh.Element{el},eDof,dV(el),full(E(el)),...
            Fem.Dim,Fem.Mesh.Node, Fem.Mesh.ShapeFnc{nn}.N,...
            Fem.Mesh.ShapeFnc{nn}.dNdxi, Fem.Mesh.ShapeFnc{nn}.W,...
            Fem.solver.sol.x, Fem.solver.sol.dx, Material.params.Rho, ...
            Material.params.Zeta, Fem.system.Gravity,...
            Material.params.Mu,  Material.params.Lambda);
    elseif isa(Material,'Yeoh') 
        nn = length(Fem.Mesh.Element{el});
        [Fe,Fb,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te,Ve,Vge,~,Re,Ue] = ...
            LocalsYHFast(Fem.Mesh.Element{el},eDof,dV(el),full(V(el)),...
            Fem.Dim,Fem.Mesh.Node, Fem.Mesh.ShapeFnc{nn}.N,...
            Fem.Mesh.ShapeFnc{nn}.dNdxi, Fem.Mesh.ShapeFnc{nn}.W,...
            Fem.solver.sol.x, Fem.solver.sol.dx, Material.params.Rho, Material.params.Zeta,...
            Fem.system.Gravity, [Material.params.C1;Material.params.C2;Material.params.C3],...
            [Material.params.D1;Material.params.D2;Material.params.D3]);   
    elseif isa(Material,'YeohIncompressible') 
            nn = length(Fem.Mesh.Element{el});
            [Fe,Fb,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te,Ve,Vge,~,Re,Ue] = ...
                LocalsYHIFast(Fem.Mesh.Element{el},eDof,dV(el),full(V(el)),...
                Fem.Dim,Fem.Mesh.Node, Fem.Mesh.ShapeFnc{nn}.N,...
                Fem.Mesh.ShapeFnc{nn}.dNdxi, Fem.Mesh.ShapeFnc{nn}.W,...
                Fem.solver.sol.x, Fem.solver.sol.dx, Material.params.Rho, Material.params.Zeta,...
                Fem.system.Gravity, [Material.params.C1;Material.params.C2;Material.params.C3]);               
    elseif isa(Material,'Mooney')
        nn = length(Fem.Element{el});
        [Fe,Fb,Qe,Me,Ce,Ke,Kte,Svme,SS,EE,Te,Ve,Vge,~,Re,Ue] = ...
            LocalsMNFast_mex(Fem.Element{el},eDof,dV(el),full(V(el)),...
            Fem.Dim,Fem.Node0,...
            Fem.ShapeFnc{nn}.N,...
            Fem.ShapeFnc{nn}.dNdxi,...
            Fem.ShapeFnc{nn}.W,...
            Fem.Utmp, Fem.dUtmp,...
            Material.Density, Material.MaterialDamping,...
            Fem.Gravity,...
            [Material.C10;Material.C01],...
            Material.K);       
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
    Fem.triplets.vj(ind3,1) = dV(el);
    Fem.triplets.l(ind3)   = Fem.Mesh.Element{el}(:);
    
    Fem.system.Potential  = Fem.system.Potential + beta*Ve;
    %Fem.PotentialG = Fem.PotentialG + beta*Vge;
    Fem.triplets.v(ind3) = Qe(:);
    
    if Fem.solver.Iteration == 1
        Fem.system.F.ElemRot{el} = Re;
        Fem.system.F.ElemStr{el} = Ue;
    end
    
    index    = index + NDof^2;
    subindex = subindex + NDof/Fem.Dim;
end

    Fem = assembleDeformationGrad(Fem);
end


Fem = assembleInternalForcesFem(Fem);
Fem.options.isAssembled = true;

end