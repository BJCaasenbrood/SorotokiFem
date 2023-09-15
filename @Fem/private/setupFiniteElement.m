function Fem = setupFiniteElement(Fem)

P = meshregularizationfilter(Fem.Mesh, Fem.topology.SpatialFilterRadius,    ...
    'Periodic',Fem.topology.Periodic);

Fem.topology.SpatialFilter = P;

if (~Fem.options.isAssembled && ~Fem.options.isNonlinear)                   ...
        || Fem.options.isNonlinear
    
    Fem.triplets.i = zeros(sum(Fem.Mesh.ElemNDof.^2),1);
    
    Fem.triplets.j  = Fem.triplets.i;
    Fem.triplets.e  = Fem.triplets.i;
    Fem.triplets.fi = Fem.triplets.i;
    Fem.triplets.k  = Fem.triplets.i;
    Fem.triplets.m  = Fem.triplets.i;
    Fem.triplets.c  = Fem.triplets.i;
    Fem.triplets.t  = Fem.triplets.i;
    Fem.triplets.fb = Fem.triplets.i;
    Fem.triplets.ft = Fem.triplets.i;
    Fem.triplets.v  = Fem.triplets.i;
    
    Fem.triplets.s  = zeros(Fem.Mesh.NNode,6);
    Fem.triplets.p  = zeros(Fem.Mesh.NNode,3);
    Fem.triplets.l  = zeros(Fem.Mesh.NNode,1);
    Fem.triplets.vj = zeros(Fem.Mesh.NNode,1);
    
end

if ~isfield(Fem.system,'ContactMesh')
    Fem = assembleContactMeshFem(Fem);
end

end