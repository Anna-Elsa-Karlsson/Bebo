function config_and_run_behavior()
% Grid entry point for behavior_analysis.m (behavioural figure + stats).
% Called as: matlab -r "config_and_run_behavior(); exit"
% No FieldTrip needed.

pipelineRoot = fileparts(mfilename('fullpath'));
addpath(pipelineRoot); % ids_ALL.m

pdat.perfFile = fullfile(pipelineRoot, '..', 'Performance_All.txt'); % BEBO2/Performance_All.txt (may hold extra IDs -- restricted to ids_ALL)
pdat.fig      = fullfile(pipelineRoot, '..', 'Plots_Behavior');

behavior_analysis(pdat);

end
