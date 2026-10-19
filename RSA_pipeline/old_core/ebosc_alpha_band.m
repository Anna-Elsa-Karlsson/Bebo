function band = ebosc_alpha_band()
% Band that defines an ALPHA RHYTHM in eBOSC (ebosc_subject.m: which
% episodes count as alpha + the range left out of the background fit;
% ebosc_results.m reads the band saved in the subject files). Kept at 8-16
% Hz -- the band the current eBOSC subject results were computed with --
% when alpha_band() (power/desync) moved to 8-20 Hz (2026-10-04). Whether
% eBOSC should use conventional alpha (8-13) and beta bands instead is an
% open decision; changing this requires rerunning the eBOSC subject step.
band = [8 16];
end
