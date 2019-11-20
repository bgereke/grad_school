classdef virmenMacro < handle
    properties
        exper;
        guiState;
        string = '';
        icon = [];
        shortcut = '';
        updateExperimentProperties = true;
        updateVariablesTable = true;
        updateWorldSketch = true;
        updateObjectProperties = true;
        updateWorldsMenu = true;
        updateWorldDrawing = true;
        updateTextureSketch = true;
        updateShapeProperties = true;
        updateTextureDrawing = true;
        trackHistory = true;
    end
    methods
        function settings(varargin)
        end
        function run(varargin)
        end
    end
end