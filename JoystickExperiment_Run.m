
clear
clear mex
% ---------- Options ----------
opts = struct();

% ---------- Control mode ----------

opts.user_controlled     = false;   % joystick drives contrast live
opts.isBinocularPlayback = false;  % if true does a prolonged set of binocular cycles
opts.enableFeedback      = true;
opts.feedbackErrorThresh = 0.25;
opts.responseLagWindowSec = [0 0.65]; % acceptable participant motor/response lag window (s)
opts.numRuns             = 2;     % cap number of runs
opts.nonius = true ;
% ---------- Control mode ----------

opts.assertTol           = 1e-10;  % float tolerance for sanity check

opts.stimulusType       = 'movie folder'; % options: 'movie folder' (default) or 'single movie'
% user-controlled mode is binocular playback only
if opts.user_controlled
    opts.isBinocularPlayback = true;
    opts.enableFeedback = false;
end

if ~islogical(opts.user_controlled) || ~isscalar(opts.user_controlled)
    error('opts.user_controlled must be a logical scalar.');
end

% ---------- Contrast mapping ----------
opts.contrast.min   = 0.0;
opts.contrast.max   = 1.0;
opts.contrast.start = 0.5;

% Which joystick control drives contrast:
% 'slider' uses state.slider01 (already normalized 0..1)
% 'x' or 'y' uses state.x/state.y (-1..1) mapped to 0..1
opts.joy.control = 'slider';
opts.joy.deadzone  = 0.00;
opts.joy.invert    = false;

% Optional shaping
opts.joy.smoothing  = 0.15;  % EMA alpha; 0=no smoothing
opts.contrast.start = 0.5;

%--------------Paths

paths = je.defaultPaths();
cd(paths.homeDir);
addpath(paths.homeDir);

% ---------- Params ----------
[stim, display] = je.joystickPsychoParams(opts);

% ---------- Session info ----------
session = je.promptSessionInfo(paths.homeDir, opts);
session.saveDir = je.ensureOutputDir(session.subjectId, paths.homeDir);

% ---------- design  timecourses ----------
if opts.isBinocularPlayback
    stim.data.runSaved = false(1, opts.numRuns);
    stim.data.created  = datestr(now);
else
    stim = je.makeJoystickPsychoConditions(stim,  opts.numRuns);
end
% ---------- PTB keyboard stuff -----------------------
KbName('UnifyKeyNames');
KbReleaseWait;

% ---------- Gamma ----------
gammaTable = je.loadGammaTable(paths.gammaTableFile);

% ---------- PTB init ----------
[display, ptb] = je.initPtb(display, stim.fix.textSizePt, gammaTable);
if opts.nonius
    [session.offsetLeft, session.offsetRight] = je.alignmentTask( ...
        'cornermatch', session.subjectId, 'eyeAdjust', 'r', ...
        'useBgPattern', 'y', 'useJoystick', 'n', ...
        'window', ptb.win, 'ifi', display.ifi, 'winRect', ptb.winRect);
end

audio = je.initFeedbackAudio(stim.tone);
cleanupObj = onCleanup(@() je.safeCleanup(ptb, audio));
% ------------ Generate stimuli------------------
stim = je.generateStimulusTimeCourses(stim, display, opts);

% ---------- Fixation geometry ----------
fixPix = je.angle2pix(display, 0.25);
stim.spatial.outerRect = CenterRect([0 0 fixPix fixPix], ptb.winRect);
stim.spatial.innerRect = CenterRect([0 0 (fixPix-10) (fixPix-10)], ptb.winRect);

% ---------- Run loop ----------
stopAll = false;
while ~stopAll
    je.abortIfEscape(ptb.keyEscape);    
    runsComplete = sum(stim.data.runSaved);
    runIndex = runsComplete + 1;
    if runIndex > numel(stim.data.runSaved)
        break;
    end

    [vidObj, movieName] = je.openMovieForRun(runIndex, session);
    stim = je.computeMovieRects( ...
        vidObj, display, stim, ptb.winRect);

    je.showInterRunScreen(ptb, session, runsComplete, numel(stim.data.runSaved), stim, display, opts);

    doNext = je.waitForNextRunOrQuit(ptb);
    if ~doNext
        stopAll = true;
        break;
    end

    [resp, timing] = je.playMovieAndCollectResponses( ...
        ptb, audio,stim, session, opts, display, vidObj, runIndex)


    fprintf('Duration: %5.2f sec\n', timing.actualSec);
    fprintf('Expected: %5.2f sec\n', timing.expectedSec);

    stim.data.response(runIndex,:) = resp;
    [~, lastPart] = fileparts(session.movieFolder);

    stim.data.blurlevel = str2double(lastPart);
    stim.data.movie{runIndex} = movieName;

    stim.data.runSaved(runIndex) = true;

    % Progressive save
    save(fullfile(session.saveDir, session.outputFileBase), "stim");
    je.abortIfEscape(ptb.keyEscape);
    if runIndex >= numel(stim.data.runSaved)
        stopAll = true;
    end
end

if opts.isBinocularPlayback
    fprintf('\nEstimated participant lag per run (binocular playback):\n');
    savedRuns = find(stim.data.runSaved);
    for k = 1:numel(savedRuns)
        runIdx = savedRuns(k);
        lagSec = je.estimateParticipantLagSec(stim, stim.data.response(runIdx,:), display, opts, runIdx);
        if isnan(lagSec)
            fprintf('  Run %d: could not estimate lag (insufficient valid samples).\n', runIdx);
        else
            fprintf('  Run %d: %0.3f s\n', runIdx, lagSec);
        end
    end
end
Screen('CloseAll');

save(fullfile(session.saveDir, session.outputFileBase), "stim", "display", "opts", "session", "ptb", "audio", "gammaTable");
