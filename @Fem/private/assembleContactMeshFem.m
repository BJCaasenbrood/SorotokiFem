function Fem = assembleContactMeshFem(Fem)

    I = (1:Fem.Mesh.NNode).';

    Fem.system.ContactMesh = struct;
    Fem.system.ContactMesh.NodeId = I;
    Fem.system.ContactMesh.qa = sort([2*I;2*I-1]);

end