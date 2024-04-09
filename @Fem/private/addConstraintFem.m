function Fem = addConstraintFem(Fem,varargin)   
for ii = 1:3:length(varargin)
    if size(varargin{ii+2},2) == 2 || size(varargin{ii+2},2) == 3 || ...
            isa(varargin{ii+2},'function_handle') || isa(varargin{ii+2},'double')

            if strcmp(varargin{ii},'Pressure')

                if ~isfield(Fem.system,'Pressure')
                    Fem.system.Pressure = [];
                end

                if ~isa(varargin{ii+2},'function_handle')
                    Fem.system.(varargin{ii}) = [Fem.system.(varargin{ii}) ; ...
                        varargin(ii+1),repmat(varargin{ii+2},[length(varargin{ii+1}),1])];
                else
                    f = varargin{ii+2};
                    Fem.system.(varargin{ii}) = [Fem.system.(varargin{ii}); {varargin{ii+1}},{f}];
                    1
                end       
            end     
        
        % if strcmpi(varargin{ii},'PressureCell')
        %     Fem.VolumetricPressure = true;
        %     if ~isa(varargin{ii+2},'function_handle')
        %         Fem.(varargin{ii}) = [Fem.(varargin{ii}); ...
        %             varargin{ii+1},repmat(varargin{ii+2},[length(varargin{ii+1}),1])];
        %     else
        %         f = varargin{ii+2};
        %         Fem.(varargin{ii}) = [Fem.(varargin{ii}); {varargin{ii+1}},{f}];
        %     end
        % elseif strcmpi(varargin{ii},'Displace')
        %     Fem.PrescribedDisplacement = true;
        %     BC = [varargin{ii+1},repmat(transpose(varargin{ii+2}(:)),...
        %         [length(varargin{ii+1}),1])];
        %     Fem.Load = [Fem.Load;BC];
        % elseif strcmpi(varargin{ii},'Contact')
        %     Fem.(varargin{ii}) = {varargin{ii+1},varargin{ii+2}};
        % elseif strcmp(varargin{ii},'Pressure')
        %     if ~isa(varargin{ii+2},'function_handle')
        %         Fem.(varargin{ii}) = [Fem.(varargin{ii}) ; ...
        %             varargin(ii+1),repmat(varargin{ii+2},[length(varargin{ii+1}),1])];
        %     else
        %         f = varargin{ii+2};
        %         Fem.(varargin{ii}) = [Fem.(varargin{ii}); {varargin{ii+1}},{f}];
        %     end
        % elseif strcmpi(varargin{ii},'Gravity')
        %     Fem.(varargin{ii}) = [varargin{ii+2}(:)];
        % else

        if strcmpi(varargin{ii},'support') || strcmpi(varargin{ii},'load') || ...
           strcmpi(varargin{ii},'tendon') || strcmpi(varargin{ii},'output') || ...
           strcmpi(varargin{ii},'spring')
            
             BC = [varargin{ii+1},repmat(transpose(varargin{ii+2}(:)),...
                 [length(varargin{ii+1}),1])];
             if isfield(Fem.system, varargin{ii})
                Fem.system.(varargin{ii}) = [Fem.system.(varargin{ii}); BC];
             else
                Fem.system.(varargin{ii}) = BC;
             end
        end
    else
        warning([varargin{ii}, ' has incorrect input'] );
    end
end
end