function Fem = addOutput(Fem,varargin)
    if isa(varargin{1},'char')
        varargin{1} = findNodeMesh(Fem.Mesh.Node,varargin{1});
    end 
    Fem = addConstraintFem(Fem,'Output',varargin{1:end});
end