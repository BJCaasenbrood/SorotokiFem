function F = assembleTendonForces(Fem)

F     = sparse(Fem.Mesh.NNode*Fem.Dim,1);   
List     = 1:Fem.Mesh.NElem;
NTendons = size(Fem.system.Tendon,1);

for jj = 1:NTendons

    % compute mean rotation matrix between elements
    id  = Fem.system.Tendon(jj,1);
    R   = Fem.system.Rotation;
    Rot = R{id};

    if Fem.Dim == 2
        Rot = Rot(1:2,1:2);
    end

    Ftend = Rot*Fem.system.Tendon(jj,2:end).';

    for ii = 1:Fem.Dim
        F(Fem.Dim*Fem.system.Tendon(jj,1)+(ii-Fem.Dim),1) = Ftend(ii);
    end
end

end