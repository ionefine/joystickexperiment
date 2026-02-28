function [offsetL, offsetR] = alignment_task(whichStimuli, sID, varargin)
%% function: alignment_task
%
% Inputs:
% - whichStimuli: a string corresponding to the stimulus you'd like to show
%
%     OPTIONS:
%     'cornermatch' : Displays two L-shapes, one in each eye, rotated such
%                     that the participant can create a cross by lining up
%                     their corners.
%     'crossinbox' : Displays a box with alignment lines in one eye, and a
%                    cross in the other eye. The task is to line up the
%                    cross inside the box.
%
% The script will ask you:

for i = 1:length(varargin)
    switch string(varargin(i))
        case 'eyeAdjust'
            eyeAdjust = string(varargin(i+1)); % can be 'l' or 'r'
        case 'useBgPattern'
            useBgPattern = string(varargin(i+1)); % can be 'y' or 'n'
        case 'addFlicker'
            addFlicker = string(varargin(i+1)); % can be 'y' or 'n'
        case 'useJoystick'
            useJoystick = string(varargin(i+1)); % can be 'y' or 'n'
    end
end

if ~exist('eyeAdjust', 'var')
    eyeAdjust = []; % ask the user
end
if ~exist('useBgPattern', 'var')
    useBgPattern = []; % ask the user
end
if ~exist('addFlicker', 'var')
    addFlicker = 'n'; % by default, no flicker - actually this isn't implemented yet so please don't use it!
end
if ~exist('useJoystick', 'var')
    useJoystick = 'n'; % by default, no joystick - use keyboard
end


%% Determine screen parameters

% obtain screen resolution
screenRes = Screen('Rect', 0);
% we are using side-by-side dual display; this means the x dimension spans
% two monitors. We have to divide the x dimension by half:
screenRes(3) = screenRes(3)/2;
% obtain the centre coordinates for the screen
screenCtr = screenRes(3:4)/2;

% Set stereomode:
stereomode = 4;
% stereomode = 4 puts the left eye image to the left, and the right eye
% image to the right (stereo mode 5 does the opposite)

%% Determine image parameters

% What stimulus will you use to run the task?
if isempty(whichStimuli)
    disp('===================== ERROR ! =====================');
    disp('= Kim changed the code, you now have to tell the  =');
    disp('= function which type of stimulus you want to use =');
    disp('=       to run this task. Check the help!         =');
    disp('===================== ERROR ! =====================');
    return;
else
    switch lower(whichStimuli)
        
        case 'cornermatch'
            % displays two L-shapes, one in each eye, rotated such that the
            % participant can create a cross by lining up their corners.
            
            % width of the lines that make up the cross
            penWidth = 10; % for Screen DrawLines (in pixels we think)
            
            % the length of the cross
            crossLength = 50; % in pixels
            
            % the size of a surrounding box presented in both eyes to aid with
            % fusion
            boxSize = 300; % in pixels
            
        case 'crossinbox'
            % Displays a box with alignment lines in one eye, and a cross in
            % the other eye. The task is to line up the cross in the box.
            
            penWidth = 10; % for Screen DrawLines (in pixels we think)
            crossLength = 50; % in pixels
            boxSize = 300; %crossLength*4; % size of the box in one eye
            
        otherwise
            
            disp('===================== ERROR ! =====================');
            disp('= I don''t know what stimulus you''re asking for. =');
            disp('=  Please check the code to determine what types  =');
            disp('=          of stimuli will be accepted.           =');
            disp('===================== ERROR ! =====================');
            return;
    end % end whichStimuli switch
end % end checking if whichStimulus is empty

% Colours
screens = Screen('Screens');
white = WhiteIndex(max(screens));
black = BlackIndex(max(screens));
grey = white / 2;

%% Keyboard/joystick

if useJoystick == 'y' % We will use the joystick
    
    % Check for the joystick and quit if it's not here
    numGamepads = Gamepad('GetNumGamepads');
    if numGamepads == 0
        disp('============================================================================')
        disp('GAMEPAD ERROR. Gamepad was not detected. Please connect gamepad and restart.')
        disp(' If error persists, restart Matlab itself, with gamepad already plugged in.');
        disp('============================================================================')
        return
    end
    
    gamepadNames = Gamepad('GetGamepadNamesFromIndices', 1:numGamepads);
    
    % select the correct index
    for i = 1:length(gamepadNames)
        if strcmpi('Thrustmaster, Inc. USB Game Controllers', gamepadNames{i})
            joystickIndex = i;
        end
    end
    joystickBuffer = 9000; % for getting near a 'zero' setting, makes joystick less sensitive
    
    joystickMiddle = 128; %
    
else % otherwise we will use the keyboard
    
    % Keyboard
    KbName('UnifyKeyNames');
    escapeKey = KbName('ESCAPE');
    leftKey = KbName('LeftArrow');
    rightKey = KbName('RightArrow');
    upKey = KbName('UpArrow');
    downKey = KbName('DownArrow');
end




%% Get input parameters from experimenter

% get participant ID

if isempty(sID)
    while 1
        sID = input('Participant ID:   ', 's');
        if ~(isempty(sID))
            break;
        end
    end
end

% get which eye to adjust
if isempty(eyeAdjust)
    while 1
        eyeAdjust = input('Adjust the (L)eft or (R)ight screen? ', 's');
        if strcmpi(eyeAdjust, 'l') || strcmpi(eyeAdjust, 'r')
            break;
        end
    end
end

if isempty(useBgPattern)
    % include a background pattern as a fusion cue?
    while 1
        useBgPattern = input('Use background pattern to aid fusion (y/n)? ', 's');
        if strcmpi(useBgPattern, 'y') || strcmpi(useBgPattern, 'n')
            break;
        end
    end
end


%% Initialize offset parameters
% Start with a zero offset
offsetL = [0 0]; % x and y offset
offsetR = [0 0];


%% Begin experiment
[window, windowRect] = Screen('OpenWindow', 0, 255/2, [], [], [], stereomode);


% Timing stuff
ifi = Screen('GetFlipInterval', window);
waitframes = 1;

%% Initialize background image, if asked for - to be developed later, KM
% % if you add images, put them into "bgimages" folder
% % make sure they are
% if useBgImg == 1
%
%     sigma = .5;
%     imgs = dir(['bgimages' filesep]); % get list of images
%     imgs = imgs(~ismember({imgs.name},{'.','..'})); % remove . and ..
%     for i = 1:length(imgs)
%         tempimg = imread([imgs(i).folder filesep imgs(i).name]);
%         [x,y] = meshgrid(linspace(-1, 1, min(size(tempimg(:,:,1)))));
%         for j = 1:size(tempimg,3)
%
%             tempgauss = padarray(exp(-(x.^2+y.^2)/sigma^2), (max(size(tempimg(:,:,1)))-min(size(tempimg(:,:,1))))/2, 0, 'both');
%             tempgauss = -tempgauss+1;
%             if size(tempgauss,1) ~= size(tempimg, 1)
%                 tempgauss = tempgauss';
%             end
%             tempfilter(:,:, j) = (tempgauss).*double(tempimg(:,:,j));
%         end
%
%         imgTex(i) = Screen('MakeTexture', window, tempimg);
%     end
%
% end %end useBgImg


%% Initialize bacground pattern, if asked for
nPatternDots = 800;
maxWindow = screenRes(3);
patternDots = randi(maxWindow, 2, nPatternDots);
patternDots = patternDots - maxWindow/2;
dotSize = 30;
dotCol = grey*.5;

%% Flicker parameters if asked for
if addFlicker == 'y'
    flickDurMsec = 5000;
    numSecs = 1; % number of seconds for each frame on 'flip'
    numFrames = round(numSecs / ifi);
    waitframes = numFrames;
end

%% Run program
vbl = Screen('flip', window);

HideCursor;
endTask = 0;
ListenChar(2);
tic;
while endTask == 0
    
    %% Select left-eye image buffer for drawing:
    Screen('SelectStereoDrawBuffer', window, 0); % selecting stereobuffer 0 fills the left eye (for stereomode 4)
    
    % If you are using a background stimulus:
    if useBgPattern == 'y'
        Screen('DrawDots', window, patternDots, dotSize, dotCol, screenCtr+offsetL, 0);
    end
    
    % Fill screen with grey
    Screen('FillRect', window, grey, [screenCtr(1)-boxSize/2+offsetL(1), screenCtr(2)-boxSize/2+offsetL(2), screenCtr(1)+boxSize/2+offsetL(1), screenCtr(2)+boxSize/2+offsetL(2)], penWidth);
    
    % Fill fixation box
    Screen('FrameRect', window, black, [screenCtr(1)-boxSize/2+offsetL(1), screenCtr(2)-boxSize/2+offsetL(2), screenCtr(1)+boxSize/2+offsetL(1), screenCtr(2)+boxSize/2+offsetL(2)], penWidth);
    
    % Draw the fixation:
    
    switch lower(whichStimuli)
        
        case 'cornermatch'
            % Draw left stim:
            Screen('DrawLines', window, [-crossLength, 0, 0, 0; 0, 0, 0, crossLength], penWidth, black, screenCtr+offsetL) % nonius cross
            Screen('DrawDots', window, [0 0], penWidth*1.5, black, screenCtr+offsetL, 1);
            Screen('DrawDots', window, [0 0], penWidth*.7, [255 0 0], screenCtr+offsetL, 1);

        case 'crossinbox'
            % Draw left stim:
            Screen('DrawLines', window, [-crossLength, crossLength, 0 0; 0, 0, -crossLength, crossLength], penWidth, black, screenCtr+offsetL); % cross
            Screen('DrawDots', window, [0 0], penWidth*1.5, black, screenCtr+offsetL, 1);
            Screen('DrawDots', window, [0 0], penWidth*.7, [255 0 0], screenCtr+offsetL, 1);
            mod(round(toc),2);

    end
    
    %% Select right-eye image buffer for drawing:
    Screen('SelectStereoDrawBuffer', window, 1) ; % selecting stereobuffer 1 fills the right eye (for stereomode 4)
    
    % If you are using a background stimulus:
    if useBgPattern == 'y'
        Screen('DrawDots', window, patternDots, dotSize, dotCol, screenCtr+offsetR, 0);
    end
    
    % Fill screen with grey:
    Screen('FillRect', window, grey, [screenCtr(1)-boxSize/2+offsetR(1), screenCtr(2)-boxSize/2+offsetR(2), screenCtr(1)+boxSize/2+offsetR(1), screenCtr(2)+boxSize/2+offsetR(2)], penWidth);
    
    % Fill fixation box:
    Screen('FrameRect', window, black, [screenCtr(1)-boxSize/2+offsetR(1), screenCtr(2)-boxSize/2+offsetR(2), screenCtr(1)+boxSize/2+offsetR(1), screenCtr(2)+boxSize/2+offsetR(2)], penWidth);
    
    % Draw the fixation:
    
    switch lower(whichStimuli)
        
        case 'cornermatch'
            Screen('DrawLines', window, [crossLength, 0, 0, 0; 0, 0, 0, -crossLength], penWidth, black, screenCtr+offsetR)
            Screen('DrawDots', window, [0 0], penWidth*1.5, black, screenCtr+offsetR, 1);
            Screen('DrawDots', window, [0 0], penWidth*.7, [255 0 0], screenCtr+offsetR, 1);
            
        case 'crossinbox'
            
            Screen('DrawLines', window, [-crossLength*2, -crossLength, crossLength, crossLength*2, 0, 0, 0, 0; 0, 0, 0, 0, -crossLength*2, -crossLength, crossLength, crossLength*2], penWidth, black, screenCtr+offsetR);
            Screen('DrawDots', window, [0 0], penWidth*1.5, black, screenCtr+offsetR, 1);
            Screen('DrawDots', window, [0 0], penWidth*.7, [255 0 0], screenCtr+offsetR, 1);
    end
    
    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
    
    % Get responses
    
    if useJoystick == 'y' % Check for Joystick Responses
        
        % Wait for the participant to press the button on the joystick
        while 1
            if Gamepad('GetButton', joystickIndex, 2) == 1 % Done with task - end
                disp('quit button down')
                while 1
                    if Gamepad('GetButton', joystickIndex, 2) == 0 % detect button release
                        disp('quit button released')
                        endTask = 1;
                        break
                    end
                end
                break
                
            elseif Gamepad('GetAxis', joystickIndex, 4) > (joystickMiddle + joystickBuffer)
                % move down
                [offsetL, offsetR] = move('down', eyeAdjust, offsetL, offsetR);
                break
                
            elseif  Gamepad('GetAxis', joystickIndex, 4) < (joystickMiddle-joystickBuffer)
                % move up
                [offsetL, offsetR] = move('up', eyeAdjust, offsetL, offsetR);
                break
                
            elseif  Gamepad('GetAxis', joystickIndex, 3) < (joystickMiddle -joystickBuffer)
                % move left?
                [offsetL, offsetR] = move('left', eyeAdjust, offsetL, offsetR);
                break
                
            elseif  Gamepad('GetAxis', joystickIndex, 3) > (joystickMiddle + joystickBuffer)
                % move right?
                [offsetL, offsetR] = move('right', eyeAdjust, offsetL, offsetR);
                break
            end
        end
        
    else % Otherwise, check for keyboard responses
        
        % Do keyboard stuff
        [keyIsDown, secs, keyCode, deltaSecs] = KbCheck();
        
        if keyIsDown == 1
            keyCode = find(keyCode, 1);
            
            if keyCode == leftKey
                [offsetL, offsetR] = move('left', eyeAdjust, offsetL, offsetR);
                
            elseif keyCode == rightKey
                [offsetL, offsetR] = move('right', eyeAdjust, offsetL, offsetR);
                
            elseif keyCode == upKey
                [offsetL, offsetR] = move('up', eyeAdjust, offsetL, offsetR);
                
            elseif keyCode == downKey
                [offsetL, offsetR] = move('down', eyeAdjust, offsetL, offsetR);
                
            elseif keyCode == escapeKey
                endTask = 1;
            end
            
        end
        
    end
    
end % end "endTask"!

%% Task is over - save stuff
runtime = datestr(now);
runnum = 1;
okToSave = 0;
while okToSave == 0
    if exist(['Alignment_Output' filesep sID '-' num2str(runnum) '-coords.mat'], 'file') ~= 0
        runnum = runnum + 1;
    else
        okToSave = 1;
    end
end
save (['Alignment_Output' filesep sID '-' num2str(runnum) '-coords.mat'], 'offsetL', 'offsetR', 'sID', 'runtime', 'whichStimuli');
Screen('CloseAll')
ListenChar(0);

end

function [offsetL, offsetR] = move(direction, eyeAdjust, offsetL, offsetR)

switch lower(direction)
    
    case 'up'
        if strcmpi(eyeAdjust, 'l') % if they are adjusting the left eye, add to the left eye offset
            offsetL(2) = offsetL(2) - 1;
        elseif strcmpi(eyeAdjust, 'r') % otherwise they chose right eye
            offsetR(2) = offsetR(2) - 1;
        end
        
    case 'down'
        if strcmpi(eyeAdjust, 'l') % if they are adjusting the left eye, add to the left eye offset
            offsetL(2) = offsetL(2) + 1;
        elseif strcmpi(eyeAdjust, 'r') % otherwise they chose right eye
            offsetR(2) = offsetR(2) + 1;
        end
        
    case 'left'
        if strcmpi(eyeAdjust, 'l') % if they are adjusting the left eye, add to the left eye offset
            offsetL(1) = offsetL(1) + 1;
        elseif strcmpi(eyeAdjust, 'r') % otherwise they chose right eye
            offsetR(1) = offsetR(1) + 1;
        end
        
    case 'right'
        if strcmpi(eyeAdjust, 'l') % if they are adjusting the left eye, add to the left eye offset
            offsetL(1) = offsetL(1) - 1;
        elseif strcmpi(eyeAdjust, 'r') % otherwise they chose right eye
            offsetR(1) = offsetR(1) - 1;
        end
        
end

end
