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
t.the_date = datestr(now, 'yyyymmdd'); % Grab today's date
t.the_time = datestr(now,'HHMM'); % Grab current time

% Generate unique seed for random number generator
t.my_rng_seed = sum(100*clock);
rng(t.my_rng_seed)

% Set discreen_sizeories
exp_dir =  pwd; % Set the experiment discreen_sizeory to the current path, grabbed by 'pwd'
data_dir = 'raw_data';

%%% Initialize Input %%%
% Check which device number the keyboard is assigned to
device_number = 0;
[keyboard_indices, product_names, ~] = GetKeyboardIndices;

%deviceString = 'Wired USB Keyboard'; % testing room 207
% deviceString = 'Apple Inc. Apple Keyboard'; % testing room 304
% deviceString = 'USB-HID Keyboard';  % desk keyboard
device_string = 'Apple Internal Keyboard / Trackpad'; %my laptop

for i = 1:length(product_names)
    if strcmp(product_names{i}, device_string)
        device_number = keyboard_indices(i);
        break;
    end
end

if device_number == 0
    error('No device by that name was detected');
end

% Setup key press
KbName('UnifyKeyNames');
keyPressNumbers = [KbName('LeftArrow') KbName('RightArrow')];

%% Screen Parameters

% Grab display parameters
screens = Screen('Screens'); % Grab the available screens
use_screen = max(screens); % If there are two or more displays, this will grab the most external display.
p.screen_size_pixels = Screen('rect', use_screen); % Grab the pixel resolution dimensions for the screen to be used

% Set display width
p.screen_width_cm = 28.6; % in centimeters; This is usually measured manually
% Set viewing distance
p.view_distance_cm = 57; % in centimeters; This is usually measured manually
% Calculating visual angle of chosen display
p.visual_angle_degrees = (2*atan2(p.screen_width_cm/2,p.view_distance_cm))*(180/pi);
p.pixels_per_degree = round(p.screen_size_pixels(3)/p.visual_angle_degrees); % Pixels per degree of visual angle

% Set fixation dot size
p.fixation_size_pixels = round(0.085*p.pixels_per_degree); % Calculate size of fixation in pixels

% Define any colors here
gray = 128; white = 255; black = 0;
green = [0 255 0]; red = [255 0 0];
p.fixation_color = black;

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
p.grating_size = round(4*p.pixels_per_degree);

%  Contrast
p.grating_contrast = 0.06; % initial contrast for target gratings (michelson contrast)

% Orientation
p.grating_orientation = 0;

% Spatial frequency
spatialFreq = 2;
p.grating_spatial_frequency = p.grating_size/p.pixels_per_degree * spatialFreq;

% Phase
p.grating_spatial_phase = 360;

%% WINDOW SETUP
% Open window and where to display center targets.
[window, screen_size] = Screen('OpenWindow', use_screen, gray,[],[],[],[],16);
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

%% Create Stimuli

% Create circular filter for textures
[cartesian_x_coordinates,cartesian_y_coordinates] = meshgrid((-p.grating_size/2):(p.grating_size/2)-1, (-p.grating_size/2):(p.grating_size/2)-1);
eccentricity_cartesian_coordinates = sqrt((cartesian_x_coordinates).^2+(cartesian_y_coordinates).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image

gaussian_filter = zeros(p.grating_size); 
gaussian_filter(eccentricity_cartesian_coordinates <= (p.grating_size/2)) = 1;

% Make grating texture grid
[Xt,Yt] = meshgrid(0:(p.grating_size-1), 0:(p.grating_size-1));
grating_texture_unfiltered = (sin(p.grating_spatial_frequency*2*pi/p.grating_size*(Xt.*sin(p.grating_orientation*(pi/180))+Yt.*cos(p.grating_orientation*(pi/180)))-p.grating_spatial_phase));

% Apply filter to grating texture and convert texture to rgb space

for n_contrast_level = 1:numel(p.grating_contrast)
    temp_grating_texture = (grating_texture_unfiltered .* gaussian_filter); % apply gaussian filter
    % apply contrast level and convert texture to rgb space -Luis
    % make drawable texture and store -Luis
end

%% Draw fixation
PsychHID('KbQueueCreate', device_number);
PsychHID('KbQueueStart', device_number);
escaped = 0;

% grating_coordinate_X = centerX;
% grating_coordinate_Y = centerY;
fixation_coordinate_X = centerX;
fixation_coordinate_Y = centerY;
update_coordinate = 1;
update_size = 50;

%moving stimulus vertically and horizontally
while 1
    grating_patch = CenterRectOnPoint([0 0 p.grating_size p.grating_size], fixation_coordinate_X, fixation_coordinate_Y);
    fixation_patch = CenterRectOnPoint([0 0 stim.fixationPx stim.fixationPx], fixation_coordinate_X, fixation_coordinate_Y); %defining size and location of fixation patch
    
    
    % Draw texture  -Luis
    
    % Draw Fixation
    Screen('FillOval', window, white, fixation_patch) 
    Screen('Flip', window);
    
    [pressed, firstPress] = PsychHID('KbQueueCheck', device_number);
    if pressed == 1
        whichPress = find(firstPress);
        if any(ismember(whichPress, KbName('ESCAPE')))
            Screen('CloseAll');
            error('User exited program');
            break;
        elseif any(ismember(whichPress, KbName('UpArrow')))
            fixation_coordinate_Y = fixation_coordinate_Y - update_size;
            
        elseif any(ismember(whichPress, KbName('DownArrow')))
            fixation_coordinate_Y = fixation_coordinate_Y + update_size;
        
        elseif any(ismember(whichPress, KbName('LeftArrow')))
            fixation_coordinate_X = fixation_coordinate_X - update_size;
        
        elseif any(ismember(whichPress, KbName('RightArrow')))
            fixation_coordinate_X = fixation_coordinate_X + update_size;
            
        end
    end
   
    

end

PsychHID('KbQueueStop', device_number);
PsychHID('KbQueueFlush', device_number);






