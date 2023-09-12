function F = addBasicLoadsFem(Fem)

F     = sparse(Fem.Mesh.NNode*Fem.Dim,1);   
NLoad = size(Fem.system.Load,1);

for ii = 1:NLoad

    fDof  = Fem.system.Load{ii,1};
    feval = Fem.system.Load{ii,2}(Fem.solver.Time);

    for jj = 1:Fem.Dim
        F(Fem.Dim*fDof+(jj-Fem.Dim),1) = feval(jj);
    end
end

end

