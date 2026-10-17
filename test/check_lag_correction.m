function check_lag_correction(T, ids, labelTxt)
% compare uncorrected vs corrected Isim - Csim on the 15 signal channels
post = {'P7','P5','P3','P1','Pz','P2','P4','P6','P8','PO7','PO3','POz','PO4','PO8','Oz'};
for s = ids(:)'
  L = load(fullfile(T, 'Results_v2', 'gamma_subject', sprintf('%d_gamma.mat', s)), 'rawMeanI','rawMeanC','meanI','meanC','refMeanI','refMeanC','PI','PC','label');
  ch = ismember(L.label, post);
  m = @(x) mean(reshape(x(ch,:,:),1,[]));
  fprintf('%s %d: raw I-C %.4f | ref I-C %.4f | corrected I-C %.4f | lags I %d-%d (med %d), C %d-%d (med %d)\n', labelTxt, s, ...
    m(L.rawMeanI)-m(L.rawMeanC), m(L.refMeanI)-m(L.refMeanC), m(L.meanI)-m(L.meanC), ...
    min(L.PI.lag), max(L.PI.lag), round(median(L.PI.lag)), min(L.PC.lag), max(L.PC.lag), round(median(L.PC.lag)));
end
end
