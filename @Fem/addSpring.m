function Fem = addSpring(Fem,varargin)
    if isa(varargin{1},'char')
        varargin{1} = findNodeMesh(Fem.Mesh.Node,varargin{1});
    end
    Fem = addSpringFem(Fem,varargin{1:end});
end