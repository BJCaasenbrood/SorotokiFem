function Fem = addDisplaceFem(Fem,varargin)   
    if ~isfield(Fem.system,'Displace')
        Fem.system.Load = {};
    end

    fDof = varargin{1};
    if ~isa(varargin{2},'function_handle')
        func = @(t) varargin{2} * clamp(t,0,1);
    end

    [n,~] = size(Fem.system.Load);

    Fem.system.Displace{n+1,1} = fDof;
    Fem.system.Displace{n+1,2} = func;
end