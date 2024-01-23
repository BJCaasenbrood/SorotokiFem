function Fem = addMyocyte(Fem,varargin)
    if isa(varargin{1},'char')
        varargin{1} = Fem.Mesh.findElements(varargin{1});
    elseif isempty(varargin{1})
        error('Please assign elements for dilation.');
    end   
    Fem = addMyocyteFem(Fem,varargin{1:end});
end