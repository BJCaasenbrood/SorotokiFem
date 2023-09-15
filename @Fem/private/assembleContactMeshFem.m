function Fem = assembleContactMeshFem(Fem)

    % I = unique(Fem.Mesh.geometry.Boundary);
    I = (1:Fem.Mesh.NNode).';
    % [~,I] = ismember(X,Y,'rows');

    Fem.system.ContactMesh = struct;
    Fem.system.ContactMesh.NodeId = I;
    Fem.system.ContactMesh.qa = sort([2*I;2*I-1]);

end