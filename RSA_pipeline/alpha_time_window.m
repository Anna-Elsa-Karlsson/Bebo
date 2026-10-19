function win = alpha_time_window()
% The ONE time window (s, encoding) for alpha/beta power in all analyses
% downstream of the desync test: single-trial power (glmer, alpha x pair
% type, between-person correlation), baseline-corrected comparison
% (Fig4b), gamma5Hz side analysis -- and the window of the desync cluster
% test itself (alpha_pre2post_group_level.m). Change it here only.
%
% 0.4-1.7 s (2026-10-06; before: the Ispec cluster's extent on both time
% axes, which in practice was ~0-1.7 s): the first 0.4 s are excluded
% because power at 8-25 Hz there is dominated by the stimulus-evoked
% response (wavelet smearing ~+/-0.4 s at 8 Hz), not the ongoing rhythm --
% the same reason the desync test starts at 0.4 s. This does NOT claim
% there is no alpha/beta desync before 0.4 s; it cannot be measured
% cleanly there. Channels stay the Ispec cluster channels, band =
% alpha_band(). (Not used by the cluster-split side analysis, which takes
% each cluster's own window.)
win = [0.4 1.7];
end
