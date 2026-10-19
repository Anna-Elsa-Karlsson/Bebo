function band = alpha_band()
% The ONE alpha(/beta) power band used by every power / desync analysis
% downstream of the alpha_pre2post existence test: single-trial alpha
% (pow1st/powAVG, alphaIspec -> Ispec x alpha correlation, glmer),
% baseline-corrected alpha (Fig4b, alpha figure desync in the Ispec
% cluster), cluster-split and gamma5Hz side analyses. Change it here only.
%
% 8-25 Hz (2026-10-06; was 8-16 Hz until 2026-10-04, then 8-20 Hz): covers
% the extent of the alpha_pre2post desync cluster (from 8 Hz up to ~25 Hz in
% the late part of the window) AND the eBOSC alpha + beta rhythm bands
% (8-14 + 14-25 Hz, ebosc_bands.m) in BOTH age groups (OA beta peaks reach
% ~25 Hz) -> an alpha/beta desync measure, cf. information-via-
% desynchronisation (Hanslmayr et al., 2012). Results were identical with
% 8-16 and 8-20 Hz (user), i.e. not dependent on the exact upper edge.
% On the 2 Hz wavelet grid = 8:2:24 Hz.
%
% NOT used for eBOSC rhythm detection -- eBOSC has its own band
% (ebosc_bands.m).
band = [8 25];
end
