function Fem = addGravityFem(Fem,varargin)   
    if ~isa(varargin{2},'function_handle')
        func = varargin{2};
    end
    Fem.system.Gravity = func;
end