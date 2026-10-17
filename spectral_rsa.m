function rmat = spectral_rsa(cfg, tfr1, tfr2, timeaxis)

% cfg - config-structure
% tfr1 - time-freq matrix for cond1
% tfr2 - time-freq matrix for cond2
% timeaxis - vector describing the time-axis of tfr1/2 in s
% note: if timeaxis(tfr1) != timeaxis(tfr2), give timeaxis as a cell-array
% with time-vectors for each tfr

% some general defaults
if ~isfield(cfg,'corr_type'), cfg.corr_type = 'pearson'; end % alternative: 'spearman'
if ~isfield(cfg,'time_avg'), cfg.time_avg = 0; end % averaging across time points (in s - default = 0 s)
if iscell(timeaxis), t1 = timeaxis{1}; t2 = timeaxis{2}; else t1 = timeaxis; t2 = []; end 
clear timeaxis;

[y1,x1] = size(tfr1); [y2,x2] = size(tfr1);
if ~(y1 == y2), error('both tfr data sets must have the same freq resolution \n'), end

% determine the averaging window
if cfg.time_avg
    g = t1(2)-t1(1); % determine sampling interval
    w = floor(cfg.time_avg/g); clear g;
else
    w = 0;
end
    
% % rmat = zeros(x1-w,x2-w);
% % for i1 = 1:(x1-w)
% %     
% %     v1 = [];
% %     v1(:,1) = tfr1(:,i1:i1+w); v1 = mean(v1,2); 
% %     
% %     for i2 = 1:(x2-w)
% %         
% %         v2 = [];
% %         v2(:,1) = tfr2(:,i2:i2+w); v2 = mean(v2,2);   
% %         
% %         if (sum(isnan(v1))>length(v1)-3) | (sum(isnan(v2))>length(v2)-3)
% %             rmat(i1,i2) = NaN;
% %         else
% %             r_tmp = corr(cat(2,v1,v2),'type',cfg.corr_type,'rows','complete');
% %             rmat(i1,i2) = r_tmp(1,2);
% %         end
% %     end
% % end

rmat = corr(tfr1,tfr2,'type',cfg.corr_type,'rows','pairwise');

% get Fisher-z transformed output
rmat = atanh(rmat);


