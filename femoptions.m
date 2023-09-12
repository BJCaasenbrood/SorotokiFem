classdef femoptions
    
    properties
        BdBox;
        Dimension;
        Color;
        ColorMap;
        ColorAxis;
        Display;
        LineStyle;
        isNonlinear;
        isAssembled;
        loadingFactor;
        isPrescribed;
    end
    
    methods
        function obj = femoptions
            %SDFOPTIONS Construct an instance of this class
            %   Detailed explanation goes here
            obj.Dimension       = 2;
            obj.Color           = [0.8329 0.8329 0.8329];
            obj.ColorMap        = cmap_turbo;
            obj.Display         = @plt;
            obj.ColorAxis       = [];
            obj.LineStyle       = '-';
            obj.isNonlinear     = true;
            obj.isAssembled     = false;
            obj.loadingFactor   = 1;
            obj.isPrescribed    = false;
        end
    end
end

function plt(Fem)
    % h = Fem.Mesh;
    % h.Node = Fem.Mesh.Node + meshfield(Fem, Fem.solver.sol.x);
    % h.show();
    cla;
    showVonMisesFem(Fem);
    axis(boxhull(Fem.Mesh.Node + meshfield(Fem, Fem.solver.sol.x)));
end
