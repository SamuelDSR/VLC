clear all;
close all;

addpath('utils');
%% Parameters

% Using FEC or not
FECcoding = 0;

%% Power and Noise

N0          = 10^-21;            % spectral density level of gaussian noise at the receiver
B           = 100*10^6;          % signal bandwidth B = 20MHz
totalPower  = 1;             % total available electrical power = E[x(n)_elec^2] or = E[x(k).*conj(x(k))]

SNR         = 120;               % dB, Pre Desired electrical SNR per bit before clipping at the transmitter (excluding the dc_bias)

%% DCO-OFDM params

nOfdmSymbol = 10000;               % number of OFDM symbols
FFTsize     = 128;               % fft size
cpSize      = 32;                % cyle prefix size
nSubcar     = FFTsize/2-1;       % number of subcarrier

%% LED params

led.minCurrent = 0.1;            % minimun turn-on forward current of LED
led.maxCurrent = 1;            % maximum forward current of LED
led.dcBias     = 0.5;             % DC bias

%% LED filter (frequency response of LED)
%low pass filter type impulse response h(t) = exp(-2*pi*fc*t) where fc is the 3dB cutoff frequency
ledCutoffFrequecy = 5*10^6;
ledFilterOrder = 32;             %32 is enough
ledFilterCoeff = led_lp_channel(ledFilterOrder,ledCutoffFrequecy/B,1);

%% VLC channel filter

%Dirac channel type
vlcFilterCoeff = 4*10^-6; %Dirac channel

%Responsitivity of PD (PD filter considered as a Dirac channel)
pd = 1;

%total equivalent channel
totalChannelCoeff = pd*conv(ledFilterCoeff,vlcFilterCoeff);
%totalChannelCoeff = pd*vlcFilterCoeff;

%frequency response of total equivalent channel 
totalChannelFR = abs(fft(totalChannelCoeff,FFTsize));

channelStateInformation = totalChannelFR(2:(nSubcar+1)); %Not sure, correct

%preSNR_elec = 10*log10(totalPower*vlc_filter^2/B/N0);

%% Bit and power loading (Bit and power loading algorithmes doesn't include the clipping effect)

% Waterfilling algorithms, only power loading included
[Capacity, powerAlloc] = ofdm_waterfilling(nSubcar,totalPower,channelStateInformation,B,N0);
%powerAlloc = totalPower/nSubcar*ones(nSubcar,1);
% Fixed QAM constellation, no bit loading
nBitsPerSymbol = 2;                                  % number of bit per QAM symbol 
yita = nBitsPerSymbol*nSubcar/(FFTsize+cpSize);      % Spectral efficiency
bitAlloc    = nBitsPerSymbol*ones(nSubcar,1);

% Campello-levin algorithmes

%% Simulation results pre-allocation

SNRelec     = zeros(size(SNR));        % Electrical SNR including DC-bias (At the Tx)
SNRelecNoDC = zeros(size(SNR));        % Electrical SNR excluding DC-bias (At the Tx)
SNRoptical  = zeros(size(SNR));        % Optical SNR including DC-bias    (At the Tx)
BER         = zeros(size(SNR));        % BER
PAPRdb      = zeros(size(SNR));        % PAPR, in dB

%% Simulation Routine
for k = 1:length(SNR)
    
    %% generate random signal
    
    generator = @random_bin_generator;
    inData =  vlc_source_generate(nOfdmSymbol*nBitsPerSymbol*nSubcar, generator);
    
    %% Encoding (interleaving / Forward Error Encoding)
    
    %FEC coding
    if FECcoding ~= 0
        encoder = @conv_encoder;
        [encData, CodeTrellis] = vlc_encode(inData, encoder, []);
        
        % interleaving
        tableSize = 10;
        interleaveTable = interleav_matrix(ones(1, tableSize));
        encoder = @interleaving;
        [encData, gain] = vlc_encode(encData, encoder, interleaveTable);
    else
        encData = inData;
    end
    
    %% bit and power loading QAM modulator
    
    modulator = @qam_modulator;
    qamParams{1} = bitAlloc;
    qamParams{2} = powerAlloc;
    [modData, remBits] = vlc_modulate(encData, modulator, qamParams);
    
    %% DCO-OFDM modulator
    
    dcoModParams = [nSubcar,cpSize];
    modulator = @dco_ofdm_modulator;
    [modData, blkSize] = vlc_modulate(modData,modulator,dcoModParams);
    
    %% Scaling  
    
    % P_elec = B*yita*N0*10^(SNR(k)/10);
    % signal_amp_moyen = sum(modData)/length(modData);
    % signal_var = sum(modData.^2)/length(modData);
    % signal_peak = max(modData.^2);
    % PAPRdb(k) = 10*log10(signal_peak/signal_var);
    % alpha = sqrt(P_elec/signal_var); %scaling factor for fitting to the LED modulation interval
    % modData = alpha*modData;
    
    %% LED Clipping
    
    %lamda = 1.5;
    %led.dc_bias = lamda*sqrt(P_elec)+led.min %set the dc bias according to the signal variance    
    led_clipping=@dco_ofdm_led_filter;
    txData  = vlc_led_filter(modData, led_clipping, led);
    
    %% Calculating the Average P_elec and Average P_opt at transmitter after LED filter
    
    P_elec_avg          = sum(txData.^2)/length(txData);                %signal elec after clipping including dc-bias
    P_elec_avg_nodc     = sum((txData-led.dcBias).^2)/length(txData);  %signal elec after clipping excluding dc-bias
    P_opt_avg           = sum(txData)/length(txData);
    
    SNRelec(k)          = 10*log10(P_elec_avg/(N0*B*yita));
    SNRelecNoDC(k)      = 10*log10(P_elec_avg_nodc/(N0*B*yita));
    SNRoptical(k)       = 10*log10(P_opt_avg/(N0*B*yita));
    
    %figure;
    %plot(tx_data(1:100*blk_size));
    
    %% led and vlc channel filtering
    
    rxData              = vlc_channel_filter(txData,totalChannelCoeff);
    
    %% Add gaussian noise
    
    noise_var           = B*N0;
    %detectData          = rxData + sqrt(noise_var)*randn(size(rxData));
    detectData = rxData;
    
    %% Demodulation
    
    % DCO-Demodulation
    dcoDemodParams(1)   = {totalChannelCoeff};
    dcoDemodParams(2)   = {cpSize};  
    dcoDemodParams(3)   = {nSubcar};
    demodulator         = @dco_ofdm_demodulator;
    demodData           = vlc_demodulate(detectData, demodulator, dcoDemodParams);
    
    % QAM Demodulation
    qamDemodParams(1)   = {bitAlloc};
    qamDemodParams(2)   = {powerAlloc};
    demodulator         = @qam_demodulator;
    demodData           = vlc_demodulate(demodData, demodulator, qamDemodParams);
    
    %% Decoding
    
    if FECcoding ~= 0         
        deleaving = @de_interleaving;
        deleav_data = vlc_decode(demodData, deleaving, interleaveTable);
        
        decoder = @conv_decoder;
        decoderData= vlc_decode(deleav_data, decoder, CodeTrellis);
    else
        decoderData = demodData;
    end
    
    %% System analysis   
    
    [number_of_errors,bit_error_rate] = biterr(inData(1:end-remBits),decoderData);  % remBits not transmitted, should deduted
    BER(k) = bit_error_rate;
    
end
BER

%% Plot

% figure;
% semilogy(SNRelec,BER,'r');
% 
% figure;
% semilogy(SNRelecNoDC,BER,'g');
% 
% figure;
% semilogy(SNRoptical,BER,'b');
