function Fem = addLoadSpring(Fem,varargin)   
    if ~isfield(Fem.system,'Spring')
        Fem.system.Spring = {};
    end

    fDof = varargin{1};
    if ~isa(varargin{2},'function_handle')
        func = @(t) varargin{2};
    end

    [n,~] = size(Fem.system.Spring);

    Fem.system.Spring{n+1,1} = fDof;
    Fem.system.Spring{n+1,2} = func;
end