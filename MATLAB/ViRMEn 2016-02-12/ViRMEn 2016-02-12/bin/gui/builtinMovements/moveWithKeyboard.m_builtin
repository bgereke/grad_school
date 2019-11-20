function velocity = moveWithKeyboard(vr)
% Keyboard control movement function for ViRMEn
%   Left/Right: change view angle
%   CTRL + Left/Right: move left/right
%   Up/Down: move forward/backward
%   CTRL + Up/Down: move up/down

persistent keyboardControl

if ~isfield(keyboardControl,'forward')
    keyboardControl.forward = 0;
    keyboardControl.rotation = 0;
    keyboardControl.sideways = 0;
    keyboardControl.vertical = 0;
end

linearScale = 30;
rotationScale = 2;

switch vr.keyPressed
    case 262
        if vr.modifiers == 0
            keyboardControl.rotation = -rotationScale;
        elseif vr.modifiers == 2
            keyboardControl.sideways = linearScale;
        end
    case 263
        if vr.modifiers == 0
            keyboardControl.rotation = rotationScale;
        elseif vr.modifiers == 2
            keyboardControl.sideways = -linearScale;
        end
    case 264
        if vr.modifiers == 0
            keyboardControl.forward = -linearScale;
        elseif vr.modifiers == 2
            keyboardControl.vertical = -linearScale;
        end
    case 265
        if vr.modifiers == 0
            keyboardControl.forward = linearScale;
        elseif vr.modifiers == 2
            keyboardControl.vertical = linearScale;
        end
end
switch vr.keyReleased
    case {262, 263}
        keyboardControl.rotation = 0;
        keyboardControl.sideways = 0;
    case {264, 265}
        keyboardControl.forward = 0;
        keyboardControl.vertical = 0;
end

velocity = [keyboardControl.forward*[sin(-vr.position(4)) cos(-vr.position(4))]+keyboardControl.sideways*[cos(vr.position(4)) sin(vr.position(4))] ...
    keyboardControl.vertical keyboardControl.rotation];