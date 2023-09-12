function Fem = addSupport(Fem,varargin)
    if isa(varargin{1},'char')
        varargin{1} = findNodeMesh(Fem.Mesh.Node,varargin{1});
    end
    Fem = addConstraintFem(Fem,'Support',varargin{1:end});
    [Fem.system.FreeDofs,Fem.system.Ia,Fem.system.FixedDofs] = getFreeDofFem(Fem); 
end