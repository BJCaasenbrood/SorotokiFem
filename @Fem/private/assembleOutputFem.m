function [L, kDof] = assembleOutputFem(Fem)

    L    = sparse(Fem.Mesh.NNode*Fem.Dim,1);   
    Nout = size(Fem.system.Output,1);

    for ii = 1:Nout
    
        fDof  = Fem.system.Output(ii,1);
        feval = Fem.system.Output(ii,2:end);
    
        for jj = 1:Fem.Dim
            L(Fem.Dim*fDof+(jj-Fem.Dim),1) = feval(jj);
        end
    end

    kDof = find(L ~= 0);
    
end
    