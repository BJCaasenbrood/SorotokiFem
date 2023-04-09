classdef topooptions
    
    properties
        sol;
        VolumeInfill;
        Thickness;
        Beta;
        SpatialFilterRadius;
        SpatialFilter;
        Change;
        Objective;
        Constraint;
        MaterialInterpolation;
        Periodic;
        Penal;
        MaxPenal;
        MaxChange;
        PenalStep;
        Ersatz;
    end
    
    methods
        function obj = topooptions
            %SDFOPTIONS Construct an instance of this class
            %   Detailed explanation goes here
            obj.MaterialInterpolation = 'SIMP';
            obj.VolumeInfill = 0.3;
            obj.Ersatz       = 1e-3; 
            obj.Penal        = 1;
            obj.MaxPenal     = 4;
            obj.PenalStep    = 10;
            obj.MaxChange    = 0.03;
            obj.Beta         = 1.5;
            obj.SpatialFilterRadius = 0.5;
        end
    end
end

