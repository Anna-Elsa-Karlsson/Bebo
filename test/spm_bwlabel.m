function [L, num] = spm_bwlabel(A, conn)
% LOCAL TEST STAND-IN ONLY (the local FieldTrip copy has no SPM). Uses the
% Image Processing Toolbox. The grid FieldTrip (20180709) ships the real
% spm_bwlabel; this file lives in test/ and is never on the grid path.
if ismatrix(A) && conn == 6, conn = 4; elseif ismatrix(A) && conn > 8, conn = 8; end
[L, num] = bwlabeln(A > 0, conn);
L = double(L);
end
