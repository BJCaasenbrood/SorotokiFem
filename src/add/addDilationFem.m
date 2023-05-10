function Fem = addDilationFem(Fem,varargin)   
    if ~isfield(Fem.system,'Dilation')
        Fem.system.Dilation = {};
    end

    fDof = varargin{1};
    if ~isa(varargin{2},'function_handle')
        func = @(t) varargin{2} * clamp(t,0,1);
    else
        func = varargin{2};
    end

    [n,~] = size(Fem.system.Dilation);

    Fem.system.Dilation{n+1,1} = fDof;
    Fem.system.Dilation{n+1,2} = func;
end