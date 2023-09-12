function Fem = setMaterial(Fem,varargin)

    if isa(varargin{1},'char')
        E = findNodeMesh(Fem.Mesh.Center,varargin{1});
    else
        E = varargin{1};
    end

    assert(varargin{2} <= Fem.materials.NMat, ['Assigned material Id',...
        ' exceeds the list of materials']);
    assert(max(E) <= Fem.Mesh.NElem, ['Assigned material Id',...
        ' exceeds the number of elements']);
    Fem.materials.MatElem(E) = varargin{2};
    
end