
classdef je
    % NaturalScenesJoystickExperiment
    %
    % Single-file refactor of the natural scenes joystick experiment.
    
    % Requires (external):
    %   Psychtoolbox
    %   VideoReader (MATLAB)
    %
   

    methods (Static)

        % ===================== Params =====================
        function [stim, display] = joystickPsychoParams(opts)
            % Refactor of the JoystickPsychoParams script (function-like, no workspace side effects)

            stim = struct();
            stim.temporal.DurPre = 2; % in secs
            if opts.isBinocularPlayback == true
                stim.temporal.DurPre = 2; % in secs
                stim.temporal.nBinoCycles = 2; % in cycles
                stim.temporal.DichDur = 0; % in secs: doesn't include DurPre and the BinocularPeriod
            else
                stim.temporal.nBinoCycles = 2; % in cycles
                stim.temporal.DichDur = 48; % in secs: doesn't include DurPre and the BinocularPeriod
            end
            stim.temporal.HzSlowCycle = 1/8;
            stim.temporal.HzFastCycle = 1/6;
            stim.temporal.HzBinoCycle = mean([stim.temporal.HzSlowCycle stim.temporal.HzFastCycle]);

            stim.temporal.DurSlowCycle = 1/stim.temporal.HzSlowCycle;
            stim.temporal.DurFastCycle = 1/stim.temporal.HzFastCycle;
            stim.temporal.DurBinoCycle = 1/stim.temporal.HzBinoCycle;

            stim.spatial.centerDeg   = [0 0];
            stim.spatial.sizeDeg     = [13 13];
            stim.spatial.meanContrast = 0.5;
            stim.fix.textSizePt = 40;

            stim.tone = struct();
            stim.tone.durationSec = 0.1;
            stim.tone.freqHz      = 400;
            stim.tone.sampleRate  = 14400;
            stim.tone.amplitude   = 0.15;

            display = struct();
            display.widthCm      = 70.1;
            display.distanceCm   = 136;
            display.screenIndex  = 0;
            display.stereoMode   = 4;
            display.waitFrames   = 1;

        end

        function paths = defaultPaths()
            paths = struct();
            cd ..
            startPath = pwd;   % or any default
            chosenPath = uigetdir(pwd, 'Select Experimental Folder');

            if isequal(chosenPath, 0)
                error('User cancelled folder selection.');
            end
            paths.homeDir        = chosenPath;
            paths.gammaTableFile = [chosenPath, '/LinearizedGammaTable.mat'];
            paths.noniusDir      = [chosenPath, '/Nonius'];
        end

        % ===================== Session =====================
        function session = promptSessionInfo(homeDir)
            session = struct();
            answer  = upper(inputdlg('Enter Subject ID:', ...
                'Subject Information', ...
                [1 40]));   % [rows cols]

            if isempty(answer)
                error('Experiment cancelled: no Subject ID entered.');
            end

            session.subjectId = upper(strtrim(answer{1}));

            session.loopThroughFolders = questdlg('What type of stimulus?', ...
                'Confirm', ...
                'single movie', 'movie folder', 'movie folder');

            if strcmp(session.loopThroughFolders, 'single movie')
                [fileName, folder] = uigetfile(fullfile(homeDir,'movies','*.avi'), 'Choose a movie');
                if isequal(fileName, 0)
                    error('No movie selected.');
                end
                session.movieFile   = fileName;
                session.movieFolder = folder;
            elseif strcmp(session.loopThroughFolders, 'movie folder')
                session.movieFolder = uigetdir(fullfile(pwd, 'stimuli'), 'Select movie folder');
                if isequal(session.movieFolder, 0)
                    error('User cancelled folder selection.');
                end
            else
                error('Invalid stimulus type');
            end
            session.offsetLeft = [0 0];  session.offsetRight = [0 0]; % initialize nonius values
        end

        function saveDir = ensureOutputDir(subjectId, homeDir)
            saveDir = fullfile(homeDir, 'output', subjectId);
            if ~isfolder(saveDir)
                mkdir(saveDir);
            end
        end

        % ===================== Files =====================
        function gammaTable = loadGammaTable(gammaFile)
            try
                s = load(gammaFile, 'gammaTable');
                if ~isfield(s, 'gammaTable')
                    error('Gamma file does not contain gammaTable: %s', gammaFile);
                end
                gammaTable = s.gammaTable;
            catch
                choice = questdlg( ...
                    'Calibration file not found. Using linear values instead. Continue?', ...
                    'Calibration Warning', ...
                    'Yes', 'No', 'No');

                if ~strcmp(choice, 'Yes')
                    error('Experiment cancelled by user.');
                else
                    gammaTable = repmat(linspace(0, 1, 256)', 1, 3);
                end
            end
        end


        % ===================== Movie handling =====================


        function [vidObj, movieName] = openMovieForRun(runIndex, session)

            if strcmp(session.loopThroughFolders, 'movie folder')
                files = dir([session.movieFolder, filesep, '*.avi']);
                nfiles = length(files);


                % a is the specific integer you want to match in *_a_*.avi
                a = mod(runIndex,nfiles) + 1; % change this mapping if needed

                pattern = sprintf('*_%d_*.avi', a);
                d = dir(fullfile(session.movieFolder, pattern));

                if isempty(d)
                    error('Movie file not found in "%s" matching pattern: %s', session.movieFolder, pattern);
                elseif numel(d) > 1
                    % If you want "first match wins", replace this error with: d = d(1);
                    names = strjoin({d.name}, ', ');
                    error('Multiple movie files match pattern %s in "%s": %s', pattern, session.movieFolder, names);
                end

                movieName = fullfile(d(1).folder, d(1).name);

            else
                % If not looping through folder, assume session.movieFile already specifies a file
                if isstruct(session.movieFile)
                    movieName = fullfile(session.movieFile(1).folder, session.movieFile(1).name);
                else
                    movieName = session.movieFile; % char/string path
                end
            end

            try
                vidObj = VideoReader(movieName);
            catch ME
                error('Could not open movie file "%s". VideoReader error: %s', movieName, ME.message);
            end
        end

        function stim = computeMovieRects(vidObj, display, stim, winRect)
            % natural-scenes sizing logic (aspect ratio)
            stimSizeDeg = 10 * [vidObj.Width/vidObj.Height, 1];
            framePx = je.angle2pix(display, stimSizeDeg); % expects [w h]
            stim.spatial.texRect = CenterRect([0 0 framePx], [0 0 display.resolutionPx]);
            stim.spatial.fusionRect = CenterRect([0 0 framePx], winRect);
        end

        % ===================== PTB init / cleanup =====================
        function [display, ptb] = initPtb(display, textSize, gammaTable)

            Screen('Preference', 'SkipSyncTests', 1);
            KbName('UnifyKeyNames');

            PsychImaging('PrepareConfiguration');

            % Screen selection
            if isstruct(display) && isfield(display,'screenIndex') && ~isempty(display.screenIndex)
                screenNumber = display.screenIndex;
            else
                screenNumber = max(Screen('Screens'));
            end

            % Stereo mode
            stereoMode = 0;
            if isstruct(display) && isfield(display,'stereoMode') && ~isempty(display.stereoMode)
                stereoMode = display.stereoMode;
            end

            [ptb.win, ptb.winRect] = PsychImaging('OpenWindow', screenNumber, 128, [], [], 2, stereoMode);
            display.resolutionPx = [ptb.winRect(4) ptb.winRect(3)];
            % Gamma
            if nargin >= 3 && ~isempty(gammaTable)
                Screen('LoadNormalizedGammaTable', ptb.win, gammaTable);
            end

            Screen('BlendFunction', ptb.win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            Screen('TextSize', ptb.win, textSize);

            Priority(MaxPriority(ptb.win));
            display.ifi = Screen('GetFlipInterval', ptb.win);
            display.frameRateHz = 1/display.ifi;

            ptb.flipTimestamp = Screen('Flip', ptb.win);
            HideCursor;
            ptb.keyEscape = KbName('ESCAPE');

            % -----------------------------
            % Thrustmaster / joystick input
            % -----------------------------
            ptb.input = struct();

            if ispc
                % WinJoystickMex is PTB's Windows joystick stopgap
                % (located in Psychtoolbox/PsychContributed)
                if exist('WinJoystickMex','file') ~= 3 && exist('WinJoystickMex','file') ~= 2
                    contribPath = fullfile(PsychtoolboxRoot, 'PsychContributed');
                    if exist(contribPath,'dir'), addpath(contribPath); end
                end
                if exist('WinJoystickMex','file') ~= 3 && exist('WinJoystickMex','file') ~= 2
                    error(['Windows joystick support needs WinJoystickMex (PsychContributed). ' ...
                        'Make sure Psychtoolbox is installed correctly and PsychContributed is on the path.']);
                end

                ptb.input.mode = 'winjoystickmex';
                ptb.input.joyId = 0; % often 0 = first joystick; change if needed
                % Axis mapping guess (you can change after a quick test)
                ptb.input.axisMap = struct('x',1,'y',2,'slider',3); % slider often comes in as z
            else
                % macOS/Linux: PTB Gamepad
                ptb.input.mode = 'gamepad';
                numPads = Gamepad('GetNumGamepads');
                if numPads == 0
                    error('No gamepad/joystick detected.');
                end
                names = Gamepad('GetGamepadNamesFromIndices', 1:numPads);
                ptb.joystickIndex = je.selectJoystickIndex(names);

                % Axis mapping guess (you can change after a quick test)
                ptb.input.axisMap = struct('x',1,'y',2,'slider',3);
            end
        end

        function idx = selectJoystickIndex(gamepadNames)
            idx = [];
            for i = 1:numel(gamepadNames)
                nm = lower(gamepadNames{i});
                if contains(nm, 'thrustmaster') && (contains(nm, 'pro') || contains(nm, 't.16000m') || contains(nm, '16000'))
                    idx = i;
                    break;
                end
            end
            if isempty(idx)
                error('Thrustmaster Pro joystick not found. Detected devices: %s', strjoin(gamepadNames, ', '));
            end
        end

        function audio = initFeedbackAudio(tone)
            InitializePsychSound;

            t = (1/tone.sampleRate):(1/tone.sampleRate):tone.durationSec;
            wave = tone.amplitude * sin(2*pi*tone.freqHz*t);
            wavedata = [wave; wave]; % stereo

            audio = struct();
            audio.handle = [];
            audio.toneDurationSec = tone.durationSec;
            audio.lastBeepTimeSec = 0;

            device = [];
            try
                audio.handle = PsychPortAudio('Open', device, [], 0, tone.sampleRate, 2);
            catch
                psychlasterror('reset');
                audio.handle = PsychPortAudio('Open', device, [], 0, [], 2);
            end

            PsychPortAudio('FillBuffer', audio.handle, wavedata);
        end

        function safeCleanup(ptb, audio)
            try
                if ~isempty(audio) && isfield(audio,'handle') && ~isempty(audio.handle)
                    PsychPortAudio('Close', audio.handle);
                end
            catch
            end

            try
                if ~isempty(ptb) && isfield(ptb,'win') && ptb.win > 0
                    screenClut = repmat(linspace(0,1,256)', 1, 3);
                    Screen('LoadNormalizedGammaTable', ptb.win, screenClut);
                end
            catch
            end

            try
                Screen('CloseAll');
                ShowCursor;
            catch
            end
        end

        % ===================== UI =====================
        function showInterRunScreen(ptb, session, runsComplete, totalRuns, stim, display, opts)
            if opts.user_controlled
                line1 = 'YOU control the contrast with the joystick.';
                line2 = 'Move the control to adjust contrast.';
            else
                line1 = 'The COMPUTER is controlling the contrast.';
                line2 = 'Use the Thrustmaster Pro slider to REPORT what you see.';
            end

            msg = {
                line1
                line2
                'Press ESC to quit.'
                'Press any other key to do the next run.'
                ['Runs completed: ' num2str(runsComplete) ' of ' num2str(totalRuns)]
                };

            txtx = ptb.winRect(3)/3;
            textSize = Screen('TextSize', ptb.win);
            baseY = textSize * 5;

            for i = 1:numel(msg)
                y = baseY + textSize*(2*i+3);

                Screen('SelectStereoDrawBuffer', ptb.win, 0);
                DrawFormattedText(ptb.win, msg{i}, ...
                    txtx + session.offsetLeft(1), y + session.offsetLeft(2), [], [], 1);

                Screen('SelectStereoDrawBuffer', ptb.win, 1);
                DrawFormattedText(ptb.win, msg{i}, ...
                    txtx + session.offsetRight(1), y + session.offsetRight(2), [], [], 1);
            end

            je.drawFixFrame(ptb, stim, session.offsetLeft, 0);
            je.drawFixFrame(ptb, stim, session.offsetRight, 1);

            ptb.flipTimestamp = Screen('Flip', ptb.win, ptb.flipTimestamp + ((display.waitFrames-0.5)*display.ifi));
        end

        function drawFixFrame(ptb, stim, offset, stereoBuffer)
            Screen('SelectStereoDrawBuffer', ptb.win, stereoBuffer);
            Screen('FillRect', ptb.win, 0, stim.spatial.outerRect + [offset offset]);
            Screen('FillRect', ptb.win, 255, stim.spatial.innerRect + [offset offset]);
            Screen('FrameRect', ptb.win, 0, stim.spatial.fusionRect + [offset offset], 10);
        end

        function doNext = waitForNextRunOrQuit(ptb)
            doNext = false;
            while true
                [keyIsDown, ~, keyCode] = KbCheck;
                if ~keyIsDown
                    continue;
                end
                k = find(keyCode, 1);
                if k == ptb.keyEscape
                    doNext = false;
                    KbReleaseWait;
                    return;
                else
                    doNext = true;
                    KbReleaseWait;
                    return;
                end
            end
        end

        %----------------------- joystick functions -----------------
        function x = applyDeadzone(x, dz)
            if nargin < 2 || isempty(dz), dz = 0; end
            if abs(x) < dz
                x = 0;
            else
                % rescale to keep continuity at dz
                x = sign(x) * (abs(x) - dz) / (1 - dz);
            end
        end

        function y = emaUpdate(yPrev, yNew, alpha)
            if nargin < 3 || isempty(alpha), alpha = 0; end
            if alpha <= 0
                y = yNew;
            else
                y = (1 - alpha) * yPrev + alpha * yNew;
            end
        end

        function c = joystickToContrast(state, opts)
            % Returns contrast in [opts.contrast.min, opts.contrast.max]
            ctl = opts.joy.control;

            switch lower(ctl)
                case 'slider'
                    u01 = state.slider01;  % already 0..1 in your readThrustmaster
                case 'x'
                    x = je.applyDeadzone(state.x, opts.joy.deadzone);
                    u01 = (x + 1) / 2;
                case 'y'
                    y = je.applyDeadzone(state.y, opts.joy.deadzone);
                    u01 = (y + 1) / 2;
                otherwise
                    error('Unknown opts.joy.control = %s (use slider/x/y)', ctl);
            end

            if isfield(opts.joy,'invert') && opts.joy.invert
                u01 = 1 - u01;
            end

            % Clamp
            u01 = max(0, min(1, u01));

            c = opts.contrast.min + u01 * (opts.contrast.max - opts.contrast.min);
        end

        function state = readThrustmaster(ptb)
            % Returns:
            %   state.x, state.y    in [-1, +1] (approx)
            %   state.slider01      in [0, 1] (approx; normalized)
            %   state.buttons       intentionally empty (buttons are read from keyboard)

            if strcmp(ptb.input.mode, 'winjoystickmex')
                [x, y, z, ~] = WinJoystickMex(ptb.input.joyId);

                % WinJoystickMex typically returns axes already roughly in [-1..+1],
                % but different drivers/devices can vary. This keeps it robust.
                state.x = je.clampToUnit(x);
                state.y = je.clampToUnit(y);

                % Treat z as the throttle/slider by default:
                % Convert [-1..+1] -> [0..1]
                state.slider01 = (je.clampToUnit(z) + 1) / 2;

                state.joystickButtons = buttons; % bitmask from joystick
                state.buttons = je.readKeyboardButtons();

            else
                j = ptb.joystickIndex;

                ax = ptb.input.axisMap;
                state.x = Gamepad('GetAxis', j, ax.x);
                state.y = Gamepad('GetAxis', j, ax.y);

                slider = Gamepad('GetAxis', j, ax.slider);
                state.slider01 = (je.clampToUnit(slider) + 1) / 2;

                state.joystickButtons = [];
                state.buttons = je.readKeyboardButtons();
            end
        end


        function buttons = readKeyboardButtons()
            [keyIsDown, ~, keyCode] = KbCheck;
            buttons = struct();
            buttons.any = keyIsDown;
            buttons.escape = keyCode(KbName('ESCAPE'));
            buttons.space = keyCode(KbName('space'));
            buttons.left = keyCode(KbName('LeftArrow'));
            buttons.right = keyCode(KbName('RightArrow'));
            buttons.up = keyCode(KbName('UpArrow'));
            buttons.down = keyCode(KbName('DownArrow'));
        end

        function v = clampToUnit(v)
            if isnan(v), v = 0; end
            v = max(-1, min(1, v));
        end

        % ===================== Core frame loop =====================
       function [response, timing] = playMovieAndCollectResponses( ...
        ptb, audio, stim, session, opts, display, vidObj, runIndex)

    nFrames = numel(stim.data.t);
    response = zeros(1, nFrames);   % will store user contrast (0..1)

    % Initialize joystick contrast state
    cUser = opts.contrast.start;

    tic;
    audio.lastBeepTimeSec = toc;

    for i = 1:nFrames
        if ~hasFrame(vidObj)
            vidObj.CurrentTime = 0;
            disp('Rewinding movie.');
        end

        frameRgb = readFrame(vidObj);
        img = mean(frameRgb, 3) / 128 - 1;  % preserve original normalization

        % --------- Determine contrast for this frame ----------
        if opts.user_controlled
            st = je.readThrustmaster(ptb);
            cNow = je.joystickToContrast(st, opts);
            cUser = je.emaUpdate(cUser, cNow, opts.joy.smoothing);

            % Use same contrast both eyes by default
            cL = cUser;
            cR = cUser;

            response(i) = cUser;   % store what participant set
        else
            % Computer-controlled uses your precomputed timecourses
            cL = stim.data.contrast(runIndex, i, 1);
            cR = stim.data.contrast(runIndex, i, 2);

            % If you still want a response trace in this mode, record slider:
            st = je.readThrustmaster(ptb);
            response(i) = st.slider01;
        end
        % -----------------------------------------------------

        stimL = img * cL * 128 + 128;
        stimR = img * cR * 128 + 128;

        stimL = uint8(min(max(stimL, 0), 255));
        stimR = uint8(min(max(stimR, 0), 255));

        texL = Screen('MakeTexture', ptb.win, stimL);
        texR = Screen('MakeTexture', ptb.win, stimR);

        % Left eye
        Screen('SelectStereoDrawBuffer', ptb.win, 0);
        Screen('DrawTexture', ptb.win, texL, [], stim.spatial.texRect + [session.offsetLeft session.offsetLeft]);
        Screen('FrameRect', ptb.win, 0, stim.spatial.fusionRect + [session.offsetLeft session.offsetLeft], 1);
        Screen('FillRect', ptb.win, 0, stim.spatial.outerRect + [session.offsetLeft session.offsetLeft]);
        Screen('FillRect', ptb.win, 255, stim.spatial.innerRect + [session.offsetLeft session.offsetLeft]);

        % Right eye
        Screen('SelectStereoDrawBuffer', ptb.win, 1);
        Screen('DrawTexture', ptb.win, texR, [], stim.spatial.texRect + [session.offsetRight session.offsetRight]);
        Screen('FrameRect', ptb.win, 0, stim.spatial.fusionRect + [session.offsetRight session.offsetRight], 1);
        Screen('FillRect', ptb.win, 0, stim.spatial.outerRect + [session.offsetRight session.offsetRight]);
        Screen('FillRect', ptb.win, 255, stim.spatial.innerRect + [session.offsetRight session.offsetRight]);

        Screen('Close', texL);
        Screen('Close', texR);

        if opts.enableFeedback
            audio = je.maybePlayCongruentFeedback(audio, display, response, stim, opts, runIndex, i);
        end

        ptb.flipTimestamp = Screen('Flip', ptb.win, ptb.flipTimestamp + ((display.waitFrames-0.5)*display.ifi));
        je.abortIfEscape();
    end

    timing = struct();
    timing.actualSec = toc;
    timing.expectedSec = max(stim.data.t);
end

        function audio = maybePlayCongruentFeedback(audio, display, response, stim, opts, runIndex, frameIndex)

            delayFrames = round(opts.feedbackDelay/display.ifi);
            if frameIndex <= delayFrames
                return;
            end

            delayedIndex = frameIndex - delayFrames;

            % Congruent-only gate 
            if ~stim.data.Bino_ON(runIndex, delayedIndex)
                return;
            end

            target = stim.data.contrast(runIndex, delayedIndex, 1);
            err = abs(response(frameIndex) - target);
            if err <= opts.feedbackErrorThresh
                return;
            end

            nowSec = toc;
            minGap = audio.toneDurationSec + err/2.5;
            if nowSec - audio.lastBeepTimeSec > minGap
                PsychPortAudio('Start', audio.handle);
                audio.lastBeepTimeSec = nowSec;
            end
        end

        % ===================== design timecourses =====================
        function stim = makeJoystickPsychoConditions(stim, numRuns)
            % makeJoystickPsychoConditionConfig
            %
            % Build the randomized run schedule for the joystick psycho experiment.
            %
            % Inputs
            %   hz.Slow, hz.Fast : temporal frequencies (Hz)
            %   stim.duration    : duration of dichoptic segment (seconds)
            %
            % Output
            %   config : struct describing each run:
            %       dropSlowEye : cycle index k (or 0) where SLOW modulation is dropped
            %       dropFastEye : cycle index k (or 0) where FAST modulation is dropped
            %       fastEye     : 0 = LE fast, 1 = RE fast
            %       runSaved    : bookkeeping flag
            %       created     : timestamp (string)

            % Random permutations of cycle indices, i.e. which cycle will
            % contain the dropout

            nSlowCycles = stim.temporal.DichDur/stim.temporal.DurSlowCycle;
            nFastCycles = stim.temporal.DichDur/stim.temporal.DurFastCycle;

            whichSlowPerm = randperm(nSlowCycles)';
            whichFastPerm = randperm(nFastCycles)';

            % Build base schedule:
            %   slow-drop runs followed by fast-drop runs
            tempTrialOrder = [whichSlowPerm zeros(nSlowCycles,1); ...
                zeros(nFastCycles,1) whichFastPerm];

            % Duplicate schedule to counterbalance eye assignment
            % Create an "eyes" vector to indicate which eye gets the FAST modulation.
            % Build a doubled schedule: first copy with eyes=0, second copy with eyes=1.
            %
            % fastEye == 0 -> FAST goes to left eye (LE), SLOW to right eye (RE)
            % fastEye == 1 -> FAST goes to right eye (RE), SLOW to left eye (LE)

            eyes = [zeros(size(tempTrialOrder,1),1); ...
                ones(size(tempTrialOrder,1),1)];

            % Two copies of trial order + eye assignment, then shuffle rows
            tempTrialOrder = Shuffle([repmat(tempTrialOrder,2,1) eyes], 2);

            % For each run:
            %   - dropSlowEye: 0 means no slow dropout this run; >0 indicates which slow cycle to drop
            %   - dropFastEye: 0 means no fast dropout this run; >0 indicates which fast cycle to drop

            stim.temporal.dropSlowEye = tempTrialOrder(1:numRuns,1);
            stim.temporal.dropFastEye = tempTrialOrder(1:numRuns,2);
            % 0=left is fast, 1=right is fast
            stim.temporal.fastEye     = tempTrialOrder(1:numRuns,3);

            stim.data.runSaved = false(size(tempTrialOrder(1:numRuns,1)));
            stim.data.created  = datestr(now);
        end

        function stim = generateStimulusTimeCourses(stim, display, opts)
            % generateStimulusTimeCourses
            %
            % Generate per-run, per-frame contrast timecourses and state masks.
            % RNG-free: deterministic given inputs.

            nruns = numel(stim.data.runSaved);

            % Frame counts
            nFramesPre  = round(stim.temporal.DurPre * display.frameRateHz);
            nFramesBino = round(stim.temporal.DurBinoCycle* display.frameRateHz);
            nFramesDich = round(stim.temporal.MainDur * display.frameRateHz);
            tDich = (0:(nFramesDich-1)) / display.frameRateHz;  % time vector for the dichoptic/monocular bit

            nFramesTotal = nFramesPre + nFramesBino + nFramesDich;
            tTotal =  (0:(nFramesTotal-1)) / display.frameRateHz;   % time vector including pre and bino
            stim.data.contrast = zeros(nruns, nFramesTotal, 2);
            stim.data.t = tTotal;

            % build preduration frames
            FramesPre = stim.spatial.meanContrast * ones(1, nFramesPre);

            % build binocular frames
            tBino = (0:(nFramesBino-1)) / display.frameRateHz;
            FramesBino = (sin(2*pi*stim.temporal.HzBinoCycle*tBino) + 1) / 2;

            % initialize masks, set them for pre and bino period
            stim.data.Pre_ON = false(nruns, nFramesTotal);
            if nFramesPre > 0
                stim.data.Pre_ON(:, 1:nFramesPre) = true;
            end

            stim.data.Bino_ON = false(nruns, nFramesTotal);
            if nFramesBino > 0
                stim.data.Bino_ON(:, nFramesPre+(1:nFramesBino)) = true;
            end

            stim.data.Mono_ON = false(nruns, nFramesTotal);

            % Build timecourses fior the dichoptic/monocular bit
            if ~opts.isBinocularPlayback
                for runNum = 1:nruns

                    tempFast = (sin(2*pi*tDich / stim.temporal.HzFastCycle) + 1) / 2;
                    tempSlow = (sin(2*pi*tDich / stim.temporal.HzSlowCycle) + 1) / 2;
                    Mono_ON = zeros(size(tempFast));

                    if stim.temporal.dropFastEye(runNum) > 0
                        k = stim.temporal.dropFastEye(runNum);
                        tempFast( ...
                            tDich > stim.temporal.DurFastCycle*(k-1) + 0.75*stim.temporal.DurFastCycle & ...
                            tDich < stim.temporal.DurFastCycle*k     + 0.75*stim.temporal.DurFastCycle) = 0;
                        Mono_ON( ...
                            tDich > stim.temporal.DurFastCycle*(k-1) + 0.75*stim.temporal.DurFastCycle & ...
                            tDich < stim.temporal.DurFastCycle*k     + 0.75*stim.temporal.DurFastCycle) = 1;
                    end

                    if stim.temporal.dropSlowEye(runNum) > 0
                        k = stim.temporal.dropSlowEye(runNum);
                        tempSlow( ...
                            tDich > stim.temporal.DurSlowCycle*(k-1) + 0.75*stim.temporal.DurSlowCycle & ...
                            tDich < stim.temporal.DurSlowCycle*k     + 0.75*stim.temporal.DurSlowCycle ) = 0;
                        Mono_ON( ...
                            tDich > stim.temporal.DurFastCycle*(k-1) + 0.75*stim.temporal.DurFastCycle & ...
                            tDich < stim.temporal.DurFastCycle*k     + 0.75*stim.temporal.DurFastCycle) = 1;
                    end

                    switch stim.temporal.fastEye(runNum)
                        case 0
                            LE = tempFast; RE = tempSlow;
                        case 1
                            RE = tempFast; LE = tempSlow;
                        otherwise
                            error('fastEye must be 0 or 1');
                    end

                    stim.data.contrast(runNum,:,1) = [FramesPre, FramesBino, LE];
                    stim.data.contrast(runNum,:,2) = [FramesPre, FramesBino, RE];
                    stim.data.Mono_ON(runNum, :) =  [0.*FramesPre, 0.*FramesBino, Mono_ON];
                end
                stim.data.Dich_ON = ~stim.data.Pre_ON & ~stim.data.Bino_ON & ~stim.data.Mono_ON;
            else % just binocular
                stim.data.contrast(1:nruns,:,1) = repmat([FramesPre, FramesBino], nruns, 1);
                stim.data.contrast(1:nruns,:,2) = repmat([FramesPre, FramesBino], nruns, 1);
                stim.data.Mono_ON = false(nruns, nFramesTotal);
                stim.data.Dich_ON = false(nruns, nFramesTotal);
            end
        end % end of runs

        %----------------utilities
        function pix = angle2pix(display,ang)
            %pix = angle2pix(display,ang)
            %
            %converts visual angles in degrees to pixels.
            %
            %Inputs:
            %display.dist (distance from screen (cm))
            %display.width (width of screen (cm))
            %display.resolution (number of pixels of display in horizontal direction)
            %
            %ang (visual angle)
            %
            %Warning: assumes isotropic (square) pixels

            %Written 11/1/07 gmb zre

            %Calculate pixel size
            pixSize = display.widthCm/display.resolutionPx(1);   %cm/pix
            sz = 2*display.distanceCm*tan(pi.*ang/(2*180));  %cm
            pix = round(sz/pixSize);   %pix


        end
        function abortIfEscape()
            [keyIsDown, ~, keyCode] = KbCheck;
            if keyIsDown && keyCode(KbName('ESCAPE'))
                fprintf('Escape pressed — exiting.\n');
                Screen('CloseAll');
                clear mex
                error('UserAbort:Escape', 'User aborted with Escape key');
            end
        end
    end
end
