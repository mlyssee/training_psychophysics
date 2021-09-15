%%% training_goal1.m
%%% Psychophysics basics training

%%% Written by Luis D. Ramirez & Lys Marcelin
%%% Created on August 17, 2021
%% Prepare
% Clean the workspace
clc; clear variables; close all; sca;
commandwindow; % Forces cursor to command window
Screen('Preference', 'SkipSyncTests',1); % 1 if running on my laptop.

% Grab date
t.theData = datestr(now, 'yyyymmdd'); % Grab today's date
t.timeStamp = datestr(now,'HHMM'); % Grab current time

% Generate unique seed for random number generator
t.mySeed = sum(100*clock);
rng(t.mySeed)

% Set discreen_sizeories
expDir =  pwd; % Set the experiment discreen_sizeory to the current path, grabbed by 'pwd'
dataDir = 'raw_data';

%%% Initialize Input %%%
% Check which device number the keyboard is assigned to
deviceNumber = 0;
[keyBoardIndices, productNames, ~] = GetKeyboardIndices;

%deviceString = 'Wired USB Keyboard'; % testing room 207
% deviceString = 'Apple Inc. Apple Keyboard'; % testing room 304
% deviceString = 'USB-HID Keyboard';  % desk keyboard
deviceString = 'Apple Internal Keyboard / Trackpad'; %my laptop

for i = 1:length(productNames)
    if strcmp(productNames{i}, deviceString)
        deviceNumber = keyBoardIndices(i);
        break;
    end
end

if deviceNumber == 0
    error('No device by that name was detected');
end

% Setup key press
KbName('UnifyKeyNames');
keyPressNumbers = [KbName('LeftArrow') KbName('RightArrow')];

%% Screen Parameters

% Grab display parameters
screens = Screen('Screens'); % Grab the available screens
useScreen = max(screens); % If there are two or more displays, this will grab the most external display.
p.screenWidthPix = Screen('rect', useScreen); % Grab the pixel resolution dimensions for the screen to be used

% Set display width
p.screenWidth = 28.6; % in centimeters; This is usually measured manually
% Set viewing distance
p.viewDistance = 57; % in centimeters; This is usually measured manually
% Calculating visual angle of chosen display
p.visAngle = (2*atan2(p.screenWidth/2,p.viewDistance))*(180/pi);
p.pixels_per_degree = round(p.screenWidthPix(3)/p.visAngle); % Pixels per degree of visual angle

% Set fixation dot size
p.fixation = round(0.085*p.pixels_per_degree); % Calculate size of fixation in pixels

% Define any colors here
gray = 128; white = 255; black = 0;
green = [0 255 0]; red = [255 0 0];
p.fixationColor = black;

%% Parameters
% % Fixation
% stim.fixationDeg = 0.5; % fixation size in visual degrees
% stim.fixationPx = round(stim.fixationDeg * p.pixels_per_degree); % fixation size in pixels
% 
% %Circle parameters (using radius to make grid)
% circle_radius_degrees = 7; %visual degrees
% circle_radius_pixels = round(circle_radius_degrees * p.pixels_per_degree); %pixels
% [circleX,circleY] = meshgrid(-circle_radius_pixels:circle_radius_pixels-1,-circle_radius_pixels:circle_radius_pixels-1);
% % circleX = round(circleX); circleY = round(circleY);
% circle_coordinates = round(sqrt(circleX.^2+circleY.^2));

% Grating Parameters

% Size
p.stimulusSz = round(4*p.pixels_per_degree);

%  Contrast
p.initContrast = 0.06; % initial contrast for target gratings (michelson contrast)

% Orientation
p.gratingOrientation = 0;

% Spatial frequency
spatialFreq = 2;
p.gratingSF = p.stimulusSz/p.pixels_per_degree * spatialFreq;

% Phase
p.gratingPhase = 360;

%% Create Stimuli

% make mask to create circle for the center grating
[x,y] = meshgrid((-p.stimulusSz/2):(p.stimulusSz/2)-1, (-p.stimulusSz/2):(p.stimulusSz/2)-1);
eccen = sqrt((x).^2+(y).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image
targetGaussian = zeros(p.stimulusSz); targetGaussian(eccen <= (p.stimulusSz/2)) = 1;

% Make transparency mask for aplha blending the two images
targetTransparencyMask = zeros(p.stimulusSz); targetTransparencyMask(eccen <= ((p.stimulusSz)/2)) = 255;

% make grating
[Xt,Yt] = meshgrid(0:(p.stimulusSz-1), 0:(p.stimulusSz-1));
target = (sin(p.gratingSF*2*pi/p.stimulusSz*(Xt.*sin(p.gratingOrientation*(pi/180))+Yt.*cos(p.gratingOrientation*(pi/180)))-p.gratingPhase));
% targetGrating= (target .* targetGaussian);

% Make gratings textures
targetGratings = NaN(p.stimulusSz, p.stimulusSz);
% target = (sin(p.freq*2*pi/p.stimulusSz*(Xt.*sin(p.gratingOrientation*(pi/180))+Yt.*cos(p.gratingOrientation*(pi/180)))-p.targPhase));
targetGratings(:,:) = (target .* targetGaussian);

%% WINDOW SETUP
% Open window and where to display center targets.
[window,screen_size] = Screen('OpenWindow', useScreen, gray,[],[],[],[],16);
OriginalCLUT = Screen('ReadNormalizedGammaTable', window);
%     load('linearizedCLUT.mat'); % testingR
%     Screen('LoadNormalizedGammaTable', window, linearizedCLUT); %testingR
HideCursor;
% Enable alpha blending
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Define coordinates where to draw the stimuli
centerX = screen_size(3)/2; centerY = screen_size(4)/2; % center coordinates

% Coordinates for location in deg
centerDeg = centerX/p.pixels_per_degree; % coordinates of targets in degrees

Screen('TextStyle', window, 1);
Screen('TextSize', window, 16);
t.ifi = Screen('GetFlipInterval',window); % grab screen refresh rate

% p.degPerFrame = p.degPerSec * t.ifi;
% p.degPerFrameGabor = sind(p.degTilt) .* p.degPerFrame;


%% Draw fixation
PsychHID('KbQueueCreate', deviceNumber);
PsychHID('KbQueueStart', deviceNumber);
escaped = 0;


% grating_coordinate_X = centerX;
% grating_coordinate_Y = centerY;
fixation_coordinate_X = centerX;
fixation_coordinate_Y = centerY;
update_coordinate = 1;
update_speed = 50;

%moving stimulus vertically and horizontally
while 1
      
    %--------------------%
    %       Grating Update       %
    %--------------------%
    % Make gratings textures
    targetGratings = NaN(p.targSize, p.targSize);
    target = (sin(p.freq*2*pi/p.targSize*(Xt.*sin(p.orientation*(pi/180))+Yt.*cos(p.orientation*(pi/180)))-p.targPhase));
    targetGratings(:,:) = (target .* targetGaussian);
    gratingTexture(:,:,1) = targetGratings*p.targContrast;
    gratingTexture(:,:,2) = targetTransparencyMask;
    grating_patch = CenterRectOnPoint([0 0 p.stimulusSz p.stimulusSz], fixation_coordinate_X, fixation_coordinate_Y);
    
%     fixation_patch = CenterRectOnPoint([0 0 stim.fixationPx stim.fixationPx], fixation_coordinate_X, fixation_coordinate_Y); %defining size and location of stim patch
    
    % Draw Fixation
    Screen('FillOval', window, white, grating_patch)
    Screen('Flip', window);
    
    [pressed, firstPress] = PsychHID('KbQueueCheck', deviceNumber);
    if pressed == 1
        whichPress = find(firstPress);
        if any(ismember(whichPress, KbName('ESCAPE')))
            Screen('CloseAll');
            error('User exited program');
            break;
        elseif any(ismember(whichPress, KbName('UpArrow')))
            fixation_coordinate_Y = fixation_coordinate_Y - update_speed;
            
        elseif any(ismember(whichPress, KbName('DownArrow')))
            fixation_coordinate_Y = fixation_coordinate_Y + update_speed;
        
        elseif any(ismember(whichPress, KbName('LeftArrow')))
            fixation_coordinate_X = fixation_coordinate_X - update_speed;
        
        elseif any(ismember(whichPress, KbName('RightArrow')))
            fixation_coordinate_X = fixation_coordinate_X + update_speed;
            
        end
    end
   
    

end

PsychHID('KbQueueStop', deviceNumber);
PsychHID('KbQueueFlush', deviceNumber);






