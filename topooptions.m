classdef topooptions
    
    properties
        sol;
        Type;
        VolumeInfill;
        Thickness;
        Beta;
        SpatialFilterRadius;
        SpatialFilter;
        Iteration;
        MaxIteration;
        Change;
        Objective;
        Constraint;
        Interpolation;
        InputStiffness;
        OutputStiffness;
        Periodic;
        Repeat;
        ReflectionPlane;
        CellRepetion;
        Penal;
        MaxPenal;
        MaxChange;
        PenalStep;
        Ersatz;
        ColorMap;
    end
    
    methods
        function obj = topooptions
            %SDFOPTIONS Construct an instance of this class
            %   Detailed explanation goes here
            obj.Type = 'Compliance';
            obj.Interpolation = 'SIMP';
            obj.VolumeInfill = 0.3;
            obj.Ersatz       = 1e-3; 
            obj.Penal        = 1;
            obj.MaxPenal     = 6;
            obj.PenalStep    = 0.25;
            obj.MaxChange    = 0.2;
            obj.Beta         = 1.5;
            obj.SpatialFilterRadius = 1;
            obj.Iteration    = 1;
            obj.MaxIteration = 25;
            obj.ColorMap     = cmap_barney(-1);
            obj.InputStiffness  = 0.1;
            obj.OutputStiffness = 0.1;
        end
    end
end

