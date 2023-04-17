function Fem = addContactFem(Fem,varargin)   
    if ~isfield(Fem.system,'Contact')
        Fem.system.Load = {};
    end

    SDF = varargin{1};
    if ~isa(varargin{2},'function_handle')
        func = @(t) varargin{2} * clamp(t,0,1);
    end

    [n,~] = size(Fem.system.Load);

    Fem.system.Contact{n+1,1} = SDF;
    Fem.system.Contact{n+1,2} = zeros(Fem.Dim,1);
end