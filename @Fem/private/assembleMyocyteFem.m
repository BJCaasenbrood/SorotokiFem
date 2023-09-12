function F = assembleMyocyteFem(Fem)

F     = sparse(Fem.Mesh.NNode*Fem.Dim,1);   
NLoad = size(Fem.system.Load,1);
W  = ones(Fem.NElem,1);
Ft = beta*sparse(Fem.i,1,W(Fem.e).*Fem.ft);

for ii = 1:NLoad

    fDof  = Fem.system.Myocyte{ii,1};
    feval = Fem.system.Myocyte{ii,2}(Fem.solver.Time);

    for jj = 1:Fem.Dim
        F(Fem.Dim*fDof+(jj-Fem.Dim),1) = feval(jj);
    end
end

end

