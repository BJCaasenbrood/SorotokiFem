function Fem = addGravityFem(Fem,varargin)   
    % if ~isfield(Fem.system,'Gravity')
    %     Fem.system.Gravity = [0,-9810];
    % end

    if ~isa(varargin{2},'function_handle')
        func = varargin{2};
    end

    % [n,~] = size(Fem.system.Pressure);

    Fem.system.Gravity = func;
    % Fem.system.Pressure{n+1,2} = func;
end