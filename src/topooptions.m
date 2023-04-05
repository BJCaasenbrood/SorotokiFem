classdef topooptions
    
    properties
        VolumeInfill;
        Thickness;
        SpatialFilterRadius;
        SpatialFilter;
        Change;
        Objective;
        Constraint;
        MaterialInterpolation;
    end
    
    methods
        function obj = topooptions
            %SDFOPTIONS Construct an instance of this class
            %   Detailed explanation goes here
            obj.VolumeInfill = 2;
        end
    end
end

