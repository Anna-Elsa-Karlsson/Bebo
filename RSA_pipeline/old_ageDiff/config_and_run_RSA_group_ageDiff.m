function config_and_run_RSA_group_ageDiff(contrastType, freqLoStr, freqHiStr)
% Grid entry point for the Ispec age-difference cluster test
% (rsa_group_level_ageDiff.m). Called as:
%   matlab -r "config_and_run_RSA_group_ageDiff('Ispec','30','100'); exit"
% Needs both per-age jobs (1st-level stat files); independent of fullSample.

pipelineRoot = fileparts(mfilename('fullpath'));
addpath('/home/mpib/karlsson/matlab/functions/fieldtrip-20180709/');
ft_defaults

pdat.elec = fullfile(pipelineRoot, '..', 'Elec');
pdat.stat = fullfile(pipelineRoot, '..', 'Results', 'Stats');

rsa_group_level_ageDiff(contrastType, str2double(freqLoStr), str2double(freqHiStr), pdat);

end
