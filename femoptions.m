classdef femoptions
    
    properties
        BdBox;
        Color;
        ColorMap;
        ColorAxis;
        Display;
        LineStyle;
        isNonlinear;
        isAssembled;
        loadingFactor;
        isPrescribed;
        VoidTolerance;
        isForceContactDamping;
    end
    
    methods
        function obj = femoptions
            %SDFOPTIONS Construct an instance of this class
            %   Detailed explanation goes here
            obj.Color           = mean(cmap_barney,1);
            obj.ColorMap        = cmap_turbo;
            obj.Display         = @plt;
            obj.ColorAxis       = [];
            obj.LineStyle       = '-';
            obj.isNonlinear     = true;
            obj.isAssembled     = false;
            obj.loadingFactor   = 1;
            obj.isPrescribed    = false;
            obj.VoidTolerance   = 0.3;
            
            obj.isForceContactDamping = false;
        end
    end
end

function plt(Fem)
    cla;
    showVonMisesFem(Fem);
    axis(boxhull(Fem.Mesh.Node + meshfield(Fem, Fem.solver.sol.x)));
end
