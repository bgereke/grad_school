classdef virmenGuiState
    properties
        selectedWorld = 1;
        selectedObject = 0;
        selectedShape = 0;
        showTriangulation;
        triangulationColor;
        worldXLim;
        worldYLim;
        textureXLim;
        textureYLim;
        fileName;
        showWireframe;
    end
    methods
        function obj = virmenGuiState(defaults)
            obj.showTriangulation = defaults.showTriangulation;
            obj.triangulationColor = defaults.triangulationColor;
            obj.worldXLim = defaults.worldXLim;
            obj.worldYLim = defaults.worldYLim;
            obj.showWireframe = defaults.showWireframe;
        end
    end
end