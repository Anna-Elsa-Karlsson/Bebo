function s = display_smooth_sigma()
% The ONE display-smoothing width (in seconds) for every time x time
% similarity / t-map figure: rsa_plots.m (Fig3 matrices, Fig4 t-map),
% RSA_clusterSplit/cs_plots.m, RSA_gamma5Hz/g5_plots.m. 2-D Gaussian,
% DISPLAY ONLY -- statistics, values and cluster outlines are unaffected.
% 2026-10-03: 32 ms (= 4 samples at 8 ms; was 16 ms, stripes still visible).
s = 0.032;
end
