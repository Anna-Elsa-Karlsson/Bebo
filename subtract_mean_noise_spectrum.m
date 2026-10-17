function tfr_c = mwb_subtract_mean_noise_spectrum(cfg,tfr)

% tfr = fieldtrip structure derived from ft_freqanalysis (or similar) with rpt X chan X freq X time
% the function requires power estimates that will be fitted in
% log(power)-log(frequency) space
%
% regression lines are fitted with a robust regression


if ~isfield(cfg,'parameter'),  cfg.parameter  = 'powspctrm'; end
if ~isfield(cfg,'mode'),       cfg.mode       = 'avg';       end % estimate noise based in segment average (cfg.mode = 'time' --> separately for each time point)
if ~isfield(cfg,'rawpowflag'), cfg.rawpowflag = 1;           end % input is raw power and needs to be log-transfomed

% some preparations
d      = tfr.(cfg.parameter);
if cfg.rawpowflag % log-transformed power required
    d = log(d);
end
f(:,1) = tfr.freq;
tfr_c = rmfield(tfr,cfg.parameter); 
clear tfr;

[nr.trl, nr.chan,nr.freq,nr.time] = size(d);

d_c = zeros(size(d));
for i1 = 1:nr.trl
    
    for i2 = 1:nr.chan
        
        x = squeeze(d(i1,i2,:,:)); 
        
        if strcmp(cfg.mode,'avg')
            
            x_m(:,1) = squeeze(mean(x,2));
            
            % extimate the noise (x_n)
            b = robustfit(log(f),x_m);
            x_n = [ones(size(f)),log(f)]*b;
            clear b;
            
            % correct the noise
            x_n = repmat(x_n,1,nr.time);
            d_c(i1,i2,:,:) = x-x_n;
            
        end
        
    end
end

if cfg.rawpowflag % convert back to original input space
    d_c = exp(d_c);
end
tfr_c.(cfg.parameter) = d_c;

