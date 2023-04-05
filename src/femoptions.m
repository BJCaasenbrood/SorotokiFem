classdef femoptions
    
    properties
        BdBox;
        Dimension;
        Color;
        ColorMap;
        ColorAxis;
        Display;
        LineStyle;
    end
    
    methods
        function obj = femoptions
            %SDFOPTIONS Construct an instance of this class
            %   Detailed explanation goes here
            obj.Dimension = 2;
            obj.Color               = [32, 129, 191]/255;
            obj.ColorMap            = cmap_turbo;
            obj.Display             = true;
            obj.ColorMap            = cmap_turbo;
            obj.ColorAxis           = [];
            obj.LineStyle           = '-';
        end
    end
end

