% This files should be pre-executed before run CAP.mdl for initializing the 
% model parameters, using "set_param('CAP','PreLoadFcn', 'cap_model_init')
% befor run CAP model
clear all;
close all;

M = 16;  % CAP modulation order
nBitPerSymbol = log2(M); % nBitPerSymbol
UpSamplingFactor = 4; %upsampling factor for shaping filter

%% create root square raised cosine shaping filter   
FilterSpan = 8;  % Number of symbols of filter, even
FilterTaps = UpSamplingFactor*FilterSpan;

wt = 0:FilterSpan*2*pi/(FilterTaps):2*FilterSpan*pi;    % time axis for pulse generation (one period)
rolloff = 0.1;
filtDef = fdesign.pulseshaping(UpSamplingFactor,'Square Root Raised Cosine','N,Beta',FilterTaps,rolloff);
rrcFilter = design(filtDef);

% normalize the filter coefficient
%rrcFilter.Numerator = rrcFilter.Numerator/max(rrcFilter.Numerator);
%fvtool(rrcFilter);
inphase_pulse = rrcFilter.Numerator.* cos(wt);      % in-phase pulse
quadrature_pulse = rrcFilter.Numerator.* sin(wt);   % quadrature pulse
%normalize filter coefficient
inphase_pulse = inphase_pulse/max(abs(inphase_pulse));
quadrature_pulse = quadrature_pulse/max(abs(quadrature_pulse));

% calculated filter gain at receiver
inphase_inphase_pulse = conv(fliplr(inphase_pulse),inphase_pulse);
quadrature_quadrature_pulse = conv(fliplr(quadrature_pulse),quadrature_pulse);
inphase_quadrature_pulse=conv(fliplr(quadrature_pulse),inphase_pulse);
inphase_gain = inphase_inphase_pulse(FilterTaps+1);
quadrature_gain = quadrature_quadrature_pulse(FilterTaps+1);

%plot shaping filter and matched filer
plot(inphase_pulse);
hold on;
plot(quadrature_pulse);
% figure
% plot(inphase_inphase_pulse)
% hold on
% plot(quadrature_quadrature_pulse)
% plot(inphase_quadrature_pulse)

%% create Modulator and demodulator
hMod = modem.qammod(M);
hMod.InputType = 'Bit';
hMod.SymbolOrder = 'Gray';
hDemod = modem.qamdemod(hMod);

%% source generator
Nsym = 10000;
source = randi([0 1], Nsym*nBitPerSymbol,1);

%% transmitter
% Modulate the input data stream using 16-QAM
signal = modulate(hMod,source);

% Real and Complex part seperation
inphase_data = real(signal);
quadrature_data = imag(signal);

% Upsampling the signal
inphase_data_up = upsample(inphase_data,UpSamplingFactor);
quadrature_data_up = upsample(quadrature_data,UpSamplingFactor);

% Filter the data with TX Filter
inphaseFil = filter(inphase_pulse,1,inphase_data_up);
quadratureFil = filter(quadrature_pulse,1,quadrature_data_up);

% Summation
signalTx = inphaseFil-quadratureFil;

%% analyse transmitted signal property
% PAPR Evaluation
% nSymPerPaprUnit = 10;
% papr = zeros(Nsym/nSymPerPaprUnit,1);
% for ii = 1:Nsym/nSymPerPaprUnit
%     unit = signalTx((ii-1)*nSymPerPaprUnit*UpSamplingFactor+1:ii*nSymPerPaprUnit*UpSamplingFactor);
%     meanSquareValue = sum(unit.*unit)/(nSymPerPaprUnit*UpSamplingFactor);
%     peakValue = max(unit.*unit);
%     papr(ii) = peakValue/meanSquareValue;
% end
% % PAPR plot
% papr = 10*log10(papr);
% [n, x] = hist(papr,0:0.1:15);
% %plot(x,cumsum(n)/length(papr),'LineWidth',4)
% semilogy(x,1-cumsum(n)/length(papr),'LineWidth',2);
% xlabel('papr, x dB')
% ylabel('Probability, X <=x')
% title('CCDF plots of PAPR CAP-16')
% grid on
     
%% passing through channel
% assuming AWGN channel
EbNo = 30;
SNR = EbNo + 10*log10(nBitPerSymbol) - 10*log10(UpSamplingFactor);
singalNoisy = awgn(signalTx,SNR,'measured');

%% Receiver
% Filter the received signal with RX filter
inphaseRx = filter(fliplr(inphase_pulse),1,singalNoisy);
quadratureRx = filter(fliplr(quadrature_pulse),1,singalNoisy);

% Downsampling 
inphaseRx = downsample(inphaseRx,UpSamplingFactor);
quadratureRx = downsample(quadratureRx, UpSamplingFactor);

% Delay
delay = FilterSpan;
inphaseRx = inphaseRx(delay+1:end);
quadratureRx = quadratureRx(delay+1:end);
    
% Normalization
inphaseRx = inphaseRx/inphase_gain;
quadratureRx = quadratureRx/quadrature_gain;


% Summation
signalRx = inphaseRx - quadratureRx*1i;
%ploteye(signalRx,FilterSpan/2);
signalRx = demodulate(hDemod,signalRx);

[number_of_errors,bit_error_rate] = biterr(source(1:end-delay*nBitPerSymbol),signalRx)











