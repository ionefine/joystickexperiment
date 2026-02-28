function [offsetL, offsetR] = alignment_task(whichStimuli, sID, varargin)
% alignment_task
% Backwards-compatible wrapper around the refactored je.alignmentTask.
[offsetL, offsetR] = je.alignmentTask(whichStimuli, sID, varargin{:});
end
