function run_test_pipeline(base, YA, OA, opts)
% Full pipeline on simulated data (make_synthetic_data) in <base>.
% Uses the test-only spm_bwlabel stand-in in test/ (local FieldTrip has no SPM).
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), here);
if nargin < 4, opts = struct(); end
make_synthetic_data(base, YA, OA, opts);
setenv('BEBO_BASE', base); setenv('BEBO_IDS_FILE', fullfile(base, 'ids.txt'));
for s = [YA(:); OA(:)]'
    s1_gamma_subject(s);
    s2_alpha_subject(s);
end
g1_ispec_group();
g2_sme_group();
g5_gamma_power_group();
g3_extract();
g4_person_stats();
end
