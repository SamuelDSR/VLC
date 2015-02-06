clear all;
close all;
addpath('utils');
disp('-------------Setting simulation parameters---------------');
tic

%% Parameters
% Using FEC or not
FECcoding = 0;

%% Power and Noise
N0          = 10^-21;            % spectral density level of gaussian noise at the receiver
B           = 100*10^6;          % signal bandwidth B = 20MHz
totalPower  = 1e-1;              % total available electrical power = E[x(n)_elec^2] or = E[x(k).*conj(x(k))]

SNR         = 120;               % dB, Pre Desired electrical SNR per bit before clipping at the transmitter (excluding the dc_bias)

%% DCO-OFDM params
nOfdmSymbol = 1000;             % number of OFDM symbols
FFTsize     = 128;               % fft size
cpSize      = 36;                % cyle prefix size
nSubcar     = FFTsize/2-1;       % number of subcarrier

%% LED params
led.minCurrent = 0.1;            % minimun turn-on forward current of LED
led.maxCurrent = 2;              % maximum forward current of LED
led.dcBias     = 1;              % DC bias

%% LED filter (frequency response of LED)
% low pass filter type impulse response h(t) = exp(-2*pi*fc*t) where fc is the 3dB cutoff frequency
% for fc = 5MHz, fs = 100MHz, 32 taps is already enough to simulate this low pass channel (the last filter coefficient is already [[1.58915973878541e-05]])
% 16 taps of CP is enough to composent this channel because filter coefficient at 16 tap is already [[0.00242197535626728]]
ledCutoffFrequecy = 10*10^6;
ledFilterOrder = 31;             % 
ledFilterCoeff = led_lp_channel(ledFilterOrder,ledCutoffFrequecy/B,1);

plot(ledFilterCoeff);
title('Led LPF filter');
grid on;

%% VLC channel filter
%Dirac channel type
vlcFilterCoeff = 0.1*10^-5; %Dirac channel

%Responsitivity of PD (PD filter considered as a Dirac channel)
pd = 1;

%% Total equivalent channel
totalChannelCoeff = pd*conv(ledFilterCoeff,vlcFilterCoeff);
%totalChannelCoeff = ledFilterCoeff;
%frequency response of total equivalent channel
totalChannelFR = abs(fft(totalChannelCoeff,FFTsize));

channelStateInformation = totalChannelFR(2:(nSubcar+1)); %Not sure, correct

preSNR_elec = 10*log10(totalPower*vlcFilterCoeff^2/B/N0);

%% Bit and power loading (Bit and power loading algorithmes doesn't include the clipping effect)

% Waterfilling algorithms, only power loading included
[Capacity, powerAlloc] = ofdm_waterfilling(nSubcar,totalPower,channelStateInformation,B,N0);
powerAlloc = totalPower/nSubcar*ones(nSubcar,1);
nBitsPerQAM         = 4;                                              % number of bit per QAM symbol 
yita               	= nBitsPerQAM*nSubcar/(FFTsize+cpSize);           % Spectral efficiency
%bitAlloc            = nBitsPerQAM*randi([1,2],nSubcar,1);
bitAlloc            = nBitsPerQAM*ones(nSubcar,1);
bitAlloc(powerAlloc==0) = 0;                                          % don't load bits on null power subcarrier
snr                 = 10*log10(powerAlloc./((B*N0)/nSubcar));
%plot(snr);
  
% Campello-levin algorithmes
% maxOrder        = 10;
% bitAlloc        = randi([0 maxOrder],nSubcar,1);
% targetBER       = 10^-3;
% gap             = 2/3*erfcinv(targetBER)^2;
% subChannelNoise = (B*N0)/nSubcar;
% gainToNoise     = channelStateInformation.^2/subChannelNoise;
% [bitAlloc, powerAlloc] = campello_algo(bitAlloc,gap,gainToNoise,nSubcar,totalPower,maxOrder);

%% Simulation results pre-allocation
SNRelec     = zeros(size(SNR));        % Electrical SNR including DC-bias (At the Tx)
SNRelecNoDC = zeros(size(SNR));        % Electrical SNR excluding DC-bias (At the Tx)
SNRoptical  = zeros(size(SNR));        % Optical SNR including DC-bias    (At the Tx)
BER         = zeros(size(SNR));        % BER
PAPRdb      = zeros(size(SNR));        % PAPR, in dB

%% Simulation Routine
disp('-------------Simulation Routine Started---------------')

for k = 1:length(SNR)   
    %% generate random signal
    
    generator = @random_bin_generator;
    inData =  vlc_source_generate(nOfdmSymbol*sum(bitAlloc), generator);
    
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
    
%     yita                = 
%     SNRelec(k)          = 10*log10(P_elec_avg/(N0*B*yita));
%     SNRelecNoDC(k)      = 10*log10(P_elec_avg_nodc/(N0*B*yita));
%     SNRoptical(k)       = 10*log10(P_opt_avg/(N0*B*yita));
    
    %figure;
    %plot(tx_data(1:100*blk_size));
    
    %% led and vlc channel filtering   
    %rxData              = vlc_channel_filter(txData,totalChannelCoeff);
    rxData              = txData;
    
    %% Add gaussian noise 
    noise_var           = B*N0;
    %detectData          = rxData + sqrt(noise_var)*randn(size(rxData));
    detectData          = rxData;
    
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
toc

%% Plot

% figure;
% semilogy(SNRelec,BER,'r');
% 
% figure;
% semilogy(SNRelecNoDC,BER,'g');
% 
% figure;
% semilogy(SNRoptical,BER,'b');
