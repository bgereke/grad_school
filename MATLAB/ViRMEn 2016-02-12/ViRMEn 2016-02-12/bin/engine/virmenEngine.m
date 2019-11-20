function err = virmenEngine(exper)
% Virmen engine

% *************************************************************************
% Copyright 2013, Princeton University.  All rights reserved.
%
% By using this software the USER indicates that he or she has read,
% understood and will comply with the following:
%
%  --- Princeton University hereby grants USER nonexclusive permission to
% use, copy and/or modify this software for internal, noncommercial,
% research purposes only. Any distribution, including publication or
% commercial sale or license, of this software, copies of the software, its
% associated documentation and/or modifications of either is strictly
% prohibited without the prior consent of Princeton University. Title to
% copyright to this software and its associated documentation shall at all
% times remain with Princeton University.  Appropriate copyright notice
% shall be placed on all software copies, and a complete copy of this
% notice shall be included in all copies of the associated documentation.
% No right is granted to use in advertising, publicity or otherwise any
% trademark, service mark, or the name of Princeton University.
%
%  --- This software and any associated documentation is provided "as is"
%
% PRINCETON UNIVERSITY MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR
% IMPLIED, INCLUDING THOSE OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR
% PURPOSE, OR THAT  USE OF THE SOFTWARE, MODIFICATIONS, OR ASSOCIATED
% DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS OR
% OTHER INTELLECTUAL PROPERTY RIGHTS OF A THIRD PARTY.
%
% Princeton University shall not be liable under any circumstances for any
% direct, indirect, special, incidental, or consequential damages with
% respect to any claim by USER or any third party on account of or arising
% from the use, or inability to use, this software or its associated
% documentation, even if Princeton University has been advised of the
% possibility of those damages.
% *************************************************************************

% No error by default
err = -1;

% Load experiment
vr.exper = exper;
vr.code = exper.experimentCode(); %#ok<*STRNU>
[letterGrid, letterFont, letterAspectRatio] = virmenLoadFont;
[windows, transformations] = virmenLoadWindows(exper);
vr.windows = windows(1:4,:)';

% Load worlds
vr.worlds = struct([]);
for wNum = 1:length(vr.exper.worlds)
    vr.worlds{wNum} = loadVirmenWorld(vr.exper.worlds{wNum});
    if size(vr.worlds{wNum}.surface.colors,1) == 4
        vr.worlds{wNum}.surface.colors(4,isnan(vr.worlds{wNum}.surface.colors(4,:))) = 1-eps;
    end
end

% Initialize parameters
vr.experimentEnded = false;
vr.currentWorld = 1;
vr.position = vr.worlds{vr.currentWorld}.startLocation;
vr.movement = NaN(1,4);
vr.movementType = 'velocity';
vr.velocity = NaN(1,4);
vr.dt = NaN;
vr.dp = NaN(1,4);
vr.dpResolution = inf;
vr.collision = false;
vr.text = struct('string',{},'position',{},'size',{},'color',{},'window',{});
vr.plot = struct('x',{},'y',{},'color',{},'window',{});
vr.textClicked = NaN;
vr.keyPressed = NaN;
vr.keyReleased = NaN;
vr.buttonPressed = NaN;
vr.buttonReleased = NaN;
vr.modifiers = NaN;
vr.activeWindow = NaN;
vr.cursorPosition = NaN;
vr.iterations = 0;
vr.timeStarted = NaN;
vr.timeElapsed = 0;

% Initialize an OpenGL window
drawnow;
virmenOpenGLRoutines(0,windows,ismac);

% Run initialization code
try
    vr = vr.code.initialization(vr); %#ok<*NASGU>
catch ME
    drawnow;
    virmenOpenGLRoutines(2);
    err = struct;
    err.message = ME.message;
    err.stack = ME.stack(1:end-1);
    return
end

% Initialize engine
oldWorld = vr.currentWorld;
oldColorSize = size(vr.worlds{vr.currentWorld}.surface.colors,1);
drawnow;
virmenOpenGLRoutines(3,oldColorSize);
oldBackgroundColor = vr.worlds{vr.currentWorld}.backgroundColor;
drawnow;
virmenOpenGLRoutines(4,oldBackgroundColor);
vr.timeStarted = now;

% Run engine
while ~vr.experimentEnded
    % Update the number of iterations
    vr.iterations = vr.iterations + 1;
    
    % Determine time elapsed since the start of the iteration
    if vr.iterations == 1
        vr.dt = 0; % Don't move on the first time step
        firstTic = tic;
    else
        vr.dt = toc(firstTic) - vr.timeElapsed; % End timer to determine time step duration
        vr.timeElapsed = vr.timeElapsed + vr.dt;
    end
        
    % Input movement information
    try
        try
            numOutputs = nargout(vr.exper.movementFunction);
        catch %#ok<CTCH>
            numOutputs = 1;
        end
        if numOutputs == 1
            vr.movement = vr.exper.movementFunction(vr);
            vr.movementType = 'velocity';
        else
            [vr.movement, vr.movementType] = vr.exper.movementFunction(vr);
        end
    catch ME
        drawnow;
        virmenOpenGLRoutines(2);
        err = struct;
        err.message = ME.message;
        err.stack = ME.stack(1:end-1);
        return
    end
    
    % Process movement information
    switch vr.movementType(1)
        case 'v' % velocity
            vr.velocity = vr.movement;
            vr.dp = vr.velocity*vr.dt;
        case 'd' % displacement
            vr.dp = vr.movement;
            vr.velocity = vr.dp/vr.dt;
        case 'p' % position
            vr.dp = movement-vr.position;
            vr.velocity = vr.dp/vr.dt;
    end
    
    % Detect collisions with edges (continuous-time collision detection)
    if ~strcmp(vr.movementType(1),'p')
        [vr.dp(1:2), vr.collision] = virmenResolveCollisions(vr.position(1:2),vr.dp(1:2), ...
            vr.worlds{vr.currentWorld}.walls.endpoints,vr.worlds{vr.currentWorld}.walls.radius2, ...
            vr.worlds{vr.currentWorld}.walls.angle,vr.worlds{vr.currentWorld}.walls.border1, ...
            vr.worlds{vr.currentWorld}.walls.border2,vr.dpResolution);
    else
        vr.collision = false;
    end
    
    % Run custom code on each engine iteration
    try
        vr = vr.code.runtime(vr);
    catch ME
        drawnow;
        virmenOpenGLRoutines(2);
        err = struct;
        err.message = ME.message;
        err.stack = ME.stack(1:end-1);
        return
    end
    
    % Reset user input states (keyboard and mouse)
    vr.textClicked = NaN;
    vr.keyPressed = NaN;
    vr.keyReleased = NaN;
    vr.buttonPressed = NaN;
    vr.buttonReleased = NaN;
    vr.modifiers = NaN;
    vr.activeWindow = NaN;
    
    % Switch worlds, if necessary
    if vr.currentWorld ~= oldWorld
        oldWorld = vr.currentWorld;
    end
    
    % Set transparency options, if necessary
    if oldColorSize ~= size(vr.worlds{vr.currentWorld}.surface.colors,1)
        oldColorSize = size(vr.worlds{vr.currentWorld}.surface.colors,1);
        drawnow;
        virmenOpenGLRoutines(3,oldColorSize);
    end
    
    % Set world background color, if necessary
    if ~all(oldBackgroundColor==vr.worlds{vr.currentWorld}.backgroundColor)
        oldBackgroundColor = vr.worlds{vr.currentWorld}.backgroundColor;
        drawnow;
        virmenOpenGLRoutines(4,oldBackgroundColor);
    end
    
    % Update position
    vr.position = vr.position+vr.dp;
    
    % Translate+rotate coordinates and calculate distances from animal
    [vertexArray, distance] = virmenProcessCoordinates(vr.worlds{oldWorld}.surface.vertices,vr.position);
    
    % Transform 3D coordinates to 2D screen coordinates
    try
        try
            numInputs = nargin(vr.exper.transformationFunction);
        catch %#ok<CTCH>
            numInputs = 2;
        end
        if numInputs == 2
            vertexArrayTransformed = vr.exper.transformationFunction(vertexArray,vr);
        else
            vertexArrayTransformed = vr.exper.transformationFunction(vertexArray);
        end
    catch ME
        drawnow;
        virmenOpenGLRoutines(2);
        err = struct;
        err.message = ME.message;
        err.stack = ME.stack(1:end-1);
        return
    end
    
    % Number of transformations returned by the user's function
    nDim = size(vertexArrayTransformed,3);
    
    % Extract triangles visible in each transformation
    triangles = virmenVisibleTriangles(vr.worlds{oldWorld}.surface.triangulation,vertexArrayTransformed,nDim,size(vertexArrayTransformed,2),vr.worlds{oldWorld}.surface.visible);
    
    % Assign distances as the z coordinate
    for d = 1:nDim
        vertexArrayTransformed(3,:,d) = distance;
    end
    
    % Sort triangles from back to front (only when transparency is on)
    if size(vr.worlds{vr.currentWorld}.surface.colors,1)==4
        ord = virmenTrianglesDistance(distance,vr.worlds{oldWorld}.surface.triangulation);
        [~, ord] = sort(ord,'descend');
        triangles = virmenOrderTriangles(triangles,size(triangles,2),nDim,ord);
    end
    
    % Set up textboxes and plots
    if ~isempty(vr.text) || ~isempty(vr.plot)
        % Fill in texts with defaults
        for ndx = 1:length(vr.text)
            if isempty(vr.text(ndx).string)
                vr.text(ndx).string = '';
            end
            if isempty(vr.text(ndx).position)
                vr.text(ndx).position = [0 0];
            end
            if isempty(vr.text(ndx).size)
                vr.text(ndx).size = 0.03;
            end
            if isempty(vr.text(ndx).color)
                vr.text(ndx).color = [1 1 1];
            end
            if isempty(vr.text(ndx).window)
                vr.text(ndx).window = 1;
            end
        end
        
        % Fill in plots with defaults
        for ndx = 1:length(vr.plot)
            if isempty(vr.plot(ndx).x) || isempty(vr.plot(ndx).y)
                vr.plot(ndx).x = [];
                vr.plot(ndx).y = [];
            end
            if isempty(vr.plot(ndx).color)
                vr.plot(ndx).color = [1 1 1];
            end
            if isempty(vr.plot(ndx).window)
                vr.plot(ndx).window = 1;
            end
        end
    end
    
    % Render the environment
    drawnow;
    vr.cursorPosition = zeros(size(windows,2),2);
    for wind = 1:size(windows,2)
        % Determine the total number of line segments to draw
        tot = 0;
        for ndx = 1:length(vr.text)
            if vr.text(ndx).window == wind
                for s = 1:length(vr.text(ndx).string)
                    tot = tot+length(letterFont{double(vr.text(ndx).string(s))});
                end
            end
        end
        colors = zeros(6,tot);
        coords = zeros(4,tot);
        
        % Create arrays of coordinates and colors
        cnt = 0;
        for ndx = 1:length(vr.text)
            if vr.text(ndx).window == wind
                for s = 1:length(vr.text(ndx).string)
                    virmenCreateLetters(coords,colors,cnt,letterGrid,letterFont{double(vr.text(ndx).string(s))},vr.text(ndx).size,vr.text(ndx).position,s,vr.text(ndx).color);
                end
            end
        end
        
        % Attach plots to the arrays of coordinates and colors
        for ndx = 1:length(vr.plot)
            if vr.plot(ndx).window == wind
                sz = size(coords,2);
                coords(:,sz+1:sz+length(vr.plot(ndx).x)-1) = ...
                    [vr.plot(ndx).x(1:end-1); vr.plot(ndx).y(1:end-1); vr.plot(ndx).x(2:end); vr.plot(ndx).y(2:end)];
                colors([1 4],sz+1:sz+length(vr.plot(ndx).x)-1) = vr.plot(ndx).color(1);
                colors([2 5],sz+1:sz+length(vr.plot(ndx).x)-1) = vr.plot(ndx).color(2);
                colors([3 6],sz+1:sz+length(vr.plot(ndx).x)-1) = vr.plot(ndx).color(3);
            end
        end
        
        % Create an array of indices
        indices = 0:2*size(coords,2)-1;
        
        % Render the environment
        if ~isnan(transformations(wind)) && transformations(wind) <= nDim
            [keyPressed, keyReleased, buttonPressed, buttonReleased, modifiers, activeWindow, vr.cursorPosition(wind,:)] = ...
                virmenOpenGLRoutines(1,vertexArrayTransformed,triangles,vr.worlds{oldWorld}.surface.colors,coords,int32(indices),colors,wind,transformations(wind),3*size(vertexArrayTransformed,2),3*size(triangles,2));
        else
            [keyPressed, keyReleased, buttonPressed, buttonReleased, modifiers, activeWindow, vr.cursorPosition(wind,:)] = ...
                virmenOpenGLRoutines(1,[],[],[],coords,int32(indices),colors,wind,0,0,0);
        end
        
        % Process user inputs (keyboard and mouse)
        if keyPressed >= 0
            vr.keyPressed = keyPressed;
        end
        if keyReleased >= 0
            vr.keyReleased = keyReleased;
        end
        if buttonPressed >= 0
            vr.buttonPressed = buttonPressed;
        end
        if buttonReleased >= 0
            vr.buttonReleased = buttonReleased;
        end
        if modifiers >= 0
            vr.modifiers = modifiers;
        end
        if activeWindow >= 0
            vr.activeWindow = activeWindow+1;
        end
    end
    
    % Determine text position boundaries
    textBoundaries = zeros(length(vr.text),4);
    for ndx = 1:length(vr.text)
        textBoundaries(ndx,:) = [vr.text(ndx).position ...
            vr.text(ndx).position(1)+vr.text(ndx).size*length(vr.text(ndx).string) ...
            vr.text(ndx).position(2)+vr.text(ndx).size*letterAspectRatio];
    end
    
    % Determine which text box, if any, was clicked
    if ~isnan(vr.buttonPressed)
        x = 2*vr.cursorPosition(vr.activeWindow,1)/windows(4,vr.activeWindow)-windows(3,vr.activeWindow)/windows(4,vr.activeWindow);
        y = 1-2*vr.cursorPosition(vr.activeWindow,2)/windows(4,vr.activeWindow);
        for t = 1:length(vr.text)
            if vr.text(t).window == vr.activeWindow
                if x >= textBoundaries(t,1) && x <= textBoundaries(t,3) && ...
                        y >= textBoundaries(t,2) && y <= textBoundaries(t,4)
                    vr.textClicked = t;
                end
            end
        end
    end
    
    % Check if experiment if over
    if vr.keyPressed == 256 % Escape key
        vr.experimentEnded = true;
    end
end

% Display engine runtime information
disp(['Ran ' num2str(vr.iterations-1) ' iterations in ' num2str(vr.timeElapsed,4) ...
    ' s (' num2str(vr.timeElapsed*1000/(vr.iterations-1),3) ' ms/frame refresh time).']);

% Run termination code
try
    vr.code.termination(vr);
catch ME
    drawnow;
    virmenOpenGLRoutines(2);
    err = struct;
    err.message = ME.message;
    err.stack = ME.stack(1:end-1);
    return
end

% Close the window used by ViRMEn
drawnow;
virmenOpenGLRoutines(2);