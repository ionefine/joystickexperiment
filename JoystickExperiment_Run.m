
clear
% ---------- Options ----------
opts = struct();

% ---------- Control mode ----------
opts.user_controlled     = true;   % joystick drives contrast live
opts.computer_controlled = false;  % use stim.data.contrast timecourses
opts.isBinocularPlayback = true;  % if true does a prolonged set of binocular cycles
opts.enableFeedback      = true;
opts.feedbackDelaySec    = 0.5;
opts.feedbackErrorThresh = 0.2;
opts.nonius = false;
opts.feedbackDelay = .5; % assumed lag if providing feedback


opts.numRuns             = 10;     % cap number of runs
opts.assertTol           = 1e-10;  % float tolerance for sanity check

% user-controlled mode is binocular playback only
if opts.user_controlled
    opts.isBinocularPlayback = true;
end

if opts.user_controlled == opts.computer_controlled
    error('Set exactly one of opts.user_controlled or opts.computer_controlled to true.');
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
session = je.promptSessionInfo(paths.homeDir);
session.saveDir = je.ensureOutputDir(session.subjectId, paths.homeDir);

% ---------- design  timecourses ----------
if opts.isBinocularPlayback
    stim.data.runSaved = false(1, opts.numRuns);
    stim.data.created  = datestr(now);
else
    stim = je.makeJoystickPsychoConditions(stim,  opts.numRuns);
end
% ---------- Alignment ----------
if opts.nonius
    addpath(paths.noniusDir);
    [session.offsetLeft, session.offsetRight] = alignment_task( ...
        'cornermatch', session.subjectId, 'eyeAdjust', 'r', ...
        'useBgPattern', 'y', 'useJoystick', 'n');
end

% ---------- PTB keyboard stuff -----------------------
KbName('UnifyKeyNames');
ESCAPE = KbName('ESCAPE');
KbReleaseWait;
% ---------- Gamma ----------
gammaTable = je.loadGammaTable(paths.gammaTableFile);

% ---------- PTB init ----------
[display, ptb] = je.initPtb(display, stim.fix.textSizePt, gammaTable);
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
    je.abortIfEscape();    runsComplete = sum(stim.data.runSaved);
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

    % Progressive save (legacy filename pattern)
    save(fullfile(session.saveDir, [session.subjectId '_congruent_psychophysics_bandpass']), "stim");
    je.abortIfEscape();
    if runIndex >= numel(stim.data.runSaved)
        stopAll = true;
    end
end

save(fullfile(session.saveDir, [session.subjectId '_congruent_psychophysics_bandpass']), "stim", "display", "opts", "session", "ptb", "audio", "gammaTable");
