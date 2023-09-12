function Fem = assemblePartialFem(Fem)

dV  = 0*ones(Fem.Mesh.NElem,1);

if isfield(Fem.system,'Dilation')
    for kk = 1:size(Fem.system.Dilation,1)
        elDof = Fem.system.Dilation{kk,1};
        eldV  = Fem.system.Dilation{kk,2}(Fem.solver.Time);
        
        dV(elDof) = eldV;
    end
end

beta    = Fem.options.loadingFactor;
[E,~,~] = materialFieldFem(Fem);

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
            [Fe, Ve] = ...
                LocalsNHFastElastic_mex(Fem.Mesh.Element{el},eDof,dV(el),full(E(el)),...
                Fem.Dim,Fem.Mesh.Node, Fem.Mesh.ShapeFnc{nn}.N,...
                Fem.Mesh.ShapeFnc{nn}.dNdxi, Fem.Mesh.ShapeFnc{nn}.W,...
                Fem.solver.sol.x, Fem.solver.sol.dx, Material.params.Rho, ...
                Material.params.Zeta, Fem.system.Gravity,...
                Material.params.Mu,  Material.params.Lambda);
        end
        
        ind1 = index+1:index+NDof^2;               % matrix indexing
        ind2 = index+1:index+NDof;                 % vector indexing
        ind3 = subindex+1:subindex+NDof/Fem.Dim;   % dimension indexing
        
        I = repmat(eDof,1,NDof); J = I';
        Fem.triplets.e(ind1)  = el;
        Fem.triplets.i(ind1)  = I(:);
        Fem.triplets.j(ind1)  = J(:);
        Fem.triplets.fi(ind2) = Fe(:);
        
        Fem.triplets.vj(ind3,1) = dV(el);
        Fem.triplets.l(ind3) = Fem.Mesh.Element{el}(:);

        Fem.system.Potential = Fem.system.Potential + beta * Ve;
        
        % if Fem.solver.Iteration == 1
        %     Fem.system.F.ElemRot{el} = Re;
        %     Fem.system.F.ElemStr{el} = Ue;
        % end
        
        index    = index + NDof^2;
        subindex = subindex + NDof/Fem.Dim;
    end
    
    % Fem = assembleDeformationGrad(Fem);
end

qa = Fem.system.Ia;

% assemble the elastic force vector
Fem.system.fElastic  = sparse(Fem.triplets.i,1,E(Fem.triplets.e)       ...
    .*Fem.triplets.fi);

Fem.system.fElastic = Fem.system.fElastic(qa);
Fem.system.fDamping = Fem.system.Damping * Fem.solver.sol.dx(qa);

Fem.options.isAssembled = true;

end