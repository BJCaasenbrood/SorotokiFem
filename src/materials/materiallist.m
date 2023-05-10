classdef materiallist
    
    properties
        NMat;
        MatElem;
        Material;       
        ColorMap; 
    end
    
    methods
        function obj = materiallist
            %SDFOPTIONS Construct an instance of this class
            %   Detailed explanation goes here
            obj.NMat     = 0;
            obj.Material = [];
            obj.ColorMap = cpal_material;
        end
    end
end

