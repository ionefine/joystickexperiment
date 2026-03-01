Experiment that allows user to report the contrast of a stimulus presented dichoptically
Main experiment is JoystickExperiment_Run
Also requires je methods file.


## Data format conversion
Use `convert_old_data_format` to convert legacy `.mat` outputs into the current `stim.data.*` format used by `JoystickExperiment_Run` and `dynamic_contrast_allS`.

Examples:
- Convert one file in place: `convert_old_data_format('output/S01/S01_congruent_psychophysics_bandpass.mat')`
- Convert a full folder tree into a new location: `convert_old_data_format('legacy_output', 'output_converted')`
