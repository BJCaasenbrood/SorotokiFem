function [FreeDofs, Ia, FixedDofs] = getFreeDofFem(Fem)
FixedDofs = [];
NSupp     = size(Fem.system.Support,1);

if Fem.Dim == 2 && NSupp > 0
    FixedDofs = [Fem.system.Support(1:NSupp,2).*(2*Fem.system.Support(1:NSupp,1)-1);
        Fem.system.Support(1:NSupp,3).*(2*Fem.system.Support(1:NSupp,1))];

elseif Fem.Dim == 3 && NSupp > 0
    FixedDofs = [Fem.Support(1:NSupp,2).*(3*Fem.system.Support(1:NSupp,1)-2);
        Fem.system.Support(1:NSupp,3).*(3*Fem.system.Support(1:NSupp,1)-1);
        Fem.system.Support(1:NSupp,4).*(3*Fem.system.Support(1:NSupp,1))];
end

FixedDofs  = unique(FixedDofs(FixedDofs>0));
AllDofs    = 1:Fem.Dim*Fem.Mesh.NNode;
[FreeDofs] = setdiff(AllDofs,FixedDofs);
[Ia,~]     = ismember(AllDofs(:),FreeDofs(:));
end