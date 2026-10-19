function config_and_run_ebosc_results()
% Grid entry point for ebosc_results.m (per-subject measures + stats).
% Called as: matlab -r "config_and_run_ebosc_results(); exit"
pipelineRoot = fileparts(mfilename('fullpath'));
addpath(pipelineRoot);                                             % ids_ALL.m
addpath(fullfile(pipelineRoot, 'ebosc'));                          % THG_find_local_max_1D_20150117.m
pdat.subj  = fullfile(pipelineRoot, '..', 'Results_eBOSC', 'subjects');
pdat.stm   = fullfile(pipelineRoot, '..', 'Results', 'SingleTrialMerged'); % Ispec_alphaROI_config.mat (ROI channels)
pdat.ga    = fullfile(pipelineRoot, '..', 'Results', 'GA');        % per-subject Ispec
pdat.out   = fullfile(pipelineRoot, '..', 'Results_eBOSC');
ebosc_results(pdat);
end
