function [out, cap_filters] = cap_modulator(in, modulator_param)

FilterSpan = modulator_param(1);
UpSamplingFactor = modulator_param(2);
FilterTaps = UpSamplingFactor*FilterSpan;

wt = 0:FilterSpan*2*pi/(FilterTaps):2*FilterSpan*pi;    % time axis for pulse generation (one period)
rolloff = 0.1;
filtDef = fdesign.pulseshaping(UpSamplingFactor,'Square Root Raised Cosine','N,Beta',FilterTaps,rolloff);
rrcFilter = design(filtDef);
% fvtool(rrcFilter);

% normalize the filter coefficient
inphase_pulse = rrcFilter.Numerator.* cos(wt);      % in-phase pulse
quadrature_pulse = rrcFilter.Numerator.* sin(wt);   % quadrature pulse

inphase_pulse = inphase_pulse/max(abs(inphase_pulse));
quadrature_pulse = quadrature_pulse/max(abs(quadrature_pulse));

cap_filters = cell(1,3);
cap_filters(1) = mat2cell(inphase_pulse);
cap_filters(2) = mat2cell(quadrature_pulse);
cap_filters(3) = {UpSamplingFactor};

% real and imaginary separation
inphase_data = real(in);
quadrature_data = imag(in);

% Upsampling the signal
inphase_data_up = upsample(inphase_data,UpSamplingFactor);
quadrature_data_up = upsample(quadrature_data,UpSamplingFactor);

% Filter the data with TX Filter
inphase_data_up(end+1:end+FilterTaps/2) = 0; %zero padding
quadrature_data_up(end+1:end+FilterTaps/2) = 0; %zero padding

inphaseFil = filter(inphase_pulse,1,inphase_data_up);
quadratureFil = filter(quadrature_pulse,1,quadrature_data_up);

inphaseFil = inphaseFil(FilterTaps/2+1:end);   % eliminate delay
quadratureFil = quadratureFil(FilterTaps/2+1:end); %eliminate delay

% Summation
out = inphaseFil-quadratureFil;

% plot shaping filter and matched filer
% plot(inphase_pulse);
% hold on;
% plot(quadrature_pulse);

% figure
% plot(inphase_inphase_pulse)
% hold on
% plot(quadrature_quadrature_pulse)
% plot(inphase_quadrature_pulse)
