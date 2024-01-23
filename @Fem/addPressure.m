%-------------------------------------------------------- add pressure load
function Fem = addPressure(Fem,varargin)
    if isa(varargin{1},'char')
        varargin{1} = Fem.Mesh.findEdges(varargin{1});
    elseif isempty(varargin{1})
        varargin{1} = Fem.Mesh.findEdges('Hole');
    end   
    Fem = addPressureFem(Fem,varargin{1:end});
end