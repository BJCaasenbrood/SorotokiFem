function [F, Ke] = assembleSpringForcesFem(Fem)

    F     = sparse(Fem.Mesh.NNode*Fem.Dim,1);   
    Ke    = sparse(Fem.Mesh.NNode*Fem.Dim,Fem.Mesh.NNode*Fem.Dim);   
    NLoad = size(Fem.system.Spring,1);
    Utmp  = Fem.solver.sol.x;

    for ii = 1:NLoad
    
        kDof  = Fem.system.Spring{ii,1};
        keval = Fem.system.Spring{ii,2}(Fem.solver.Time);
    
        for jj = 1:Fem.Dim
            x = Utmp(Fem.Dim*kDof+(jj-Fem.Dim),1);
            F(Fem.Dim*kDof+(jj-Fem.Dim),1) = -keval(jj) * x;
            Ke(Fem.Dim*kDof+(jj-Fem.Dim), ...
               Fem.Dim*kDof+(jj-Fem.Dim)) = keval(jj);
        end
    end
    
end
    
    