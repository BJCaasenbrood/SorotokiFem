function [ge, ce] = assembleDisplaceEq(Fem)

    ge    = - Fem.solver.sol.x; 
    ce    = sparse(Fem.Mesh.NNode*Fem.Dim,1);   
    NLoad = size(Fem.system.Displace,1);
    
    for ii = 1:NLoad
    
        fDof  = Fem.system.Displace{ii,1};
        feval = Fem.system.Displace{ii,2}(Fem.solver.Time);
    
        for jj = 1:Fem.Dim
            ind = Fem.Dim*fDof+(jj-Fem.Dim);
            ge(ind,1) = ge(ind,1) + feval(jj);
            ce(ind,1) = 1.0;
        end

    end

end


