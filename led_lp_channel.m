%% generate a LED low pass channel
% n: filter order, (taps)
% fc: normalized cutoff frequency
% type
function filter = led_lp_channel(n,fc,type)
    t = 0:n;
    filter = exp(-2*pi*fc*t');
    % normalization, filter doesn't amplifier power
    filter = filter/sum(filter);
    %fvtool(filter,1);