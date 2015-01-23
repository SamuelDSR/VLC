clear all;
close all;

%% Parameters
%using FEC or not
FEC_Coding = 0;

%SNR
N0 = 10^-21;  % spectral density level of gaussian noise at the receiver
B = 20*10^6;  % signal bandwidth B = 20MHz
SNR = 100:1:120; %dB, Desired electrical SNR per bit before clipping at the transmitter (excluding the dc_bias)

%Cap params
Nsym = 200000;  % number of cap symbols

nBitPerSymbol = 4; %QAM order
yita = nBitPerSymbol; %spectral effiency

FilterSpan = 4;
UpSamplingFactor = 4;

%led params
led.dc_bias = 0.5; % A
led.min = 0.1; % A
led.max = 1; % A

%LED filter (Tx filter)
% led_filter = [];
% %VLC channel filter
% vlc_filter = [];
% channel = conv(led_filter,vlc_filter);
channel = 4*10^-6; %Dirac channel

%Rx filter
pd = 1;

%Simulation results pre-allocation
SNR_elec = zeros(size(SNR));
SNR_elec_nodc = zeros(size(SNR));
SNR_opt = zeros(size(SNR));
BER = zeros(size(SNR));

    = zeros(size(SNR));

for k = 1:length(SNR)
    
    %% generate random signal
    generator = @random_bin_generator;
    in_data =  vlc_source_generate(Nsym*nBitPerSymbol, generator);
    
    %% Encoding (interleaving / Forward Error Encoding)
    %FEC coding
    if FEC_Coding ~= 0
        encoder = @conv_encoder;
        [enc_data, CodeTrellis] = vlc_encode(in_data, encoder, []);
        
        % interleaving
        table_size = 10;
        interleave_table = interleav_matrix(ones(1, table_size));
        encoder = @interleaving;
        [enc_data, gain] = vlc_encode(enc_data, encoder, interleave_table);
    else
        enc_data = in_data;
    end
    
    %% Modulation
    %QAM modulator
    ModOrder = 2^nBitPerSymbol;
    modulator = @qam_modulator;
    [mod_data, qam_handle] = vlc_modulate(enc_data, modulator, ModOrder);
    
    %CAP modulator
    modulator_param = [FilterSpan,UpSamplingFactor];
    modulator = @cap_modulator;
    [mod_data, cap_params] = vlc_modulate(mod_data, modulator, modulator_param);
    
    %% Scaling the signal for matching the desired SNR
    P_elec = B*yita*N0*10^(SNR(k)/10);
    signal_amp_moyen = sum(mod_data)/length(mod_data);
    signal_var = sum(mod_data.^2)/length(mod_data);
    signal_peak = max(mod_data.^2);
    PAPR_dB(k) = 10*log10(signal_peak/signal_var);
    
    alpha = sqrt(P_elec/signal_var); %scaling factor for fitting to the LED modulation interval
    mod_data = alpha*mod_data;
    signal_var = sum(mod_data.^2)/length(mod_data);
    %% LED Clipping
    led_filter=@dco_ofdm_led_filter;
    tx_data  = vlc_led_filter(mod_data, led_filter,led);
    
    %% Calculating the Average P_elec and Average P_opt after LED filter
    P_elec_avg = sum(tx_data.^2)/length(tx_data);  %signal elec after clipping including dc-bias
    P_elec_avg_nodc = sum((tx_data-led.dc_bias).^2)/length(tx_data); %signal elec after clipping excluding dc-bias
    P_opt_avg = sum(tx_data)/length(tx_data);
    
    SNR_elec(k) = 10*log10(P_elec_avg/(N0*B*yita));
    SNR_elec_nodc(k) = 10*log10(P_elec_avg_nodc/(N0*B*yita));
    SNR_opt(k) = 10*log10(P_opt_avg/(N0*B*yita));
    
    %% channel filtering
    rx_data  = vlc_channel_filter(tx_data, channel);
    
    %% detection using receiver filter
    detect_data = vlc_pd_filter(rx_data, pd);
    
    %% Add gaussian noise
    noise_var = B*N0;
    detect_data = detect_data + sqrt(noise_var)*randn(size(detect_data));
    
    %% Demodulation
    cap_params(4) = {[alpha,channel*pd]};  %attenuation, normally acknowledged at Rx by estimation of channel (Ex: using pilot)
    cap_params(5) = {led.dc_bias};         %dc-bias added at Tx, normally removed at Rx using a High Pass filter
    
    demodulator = @cap_demodulator;
    demod_data  = vlc_demodulate(detect_data, demodulator, cap_params);
    
    demodulator = @qam_demodulator;
    demod_data  = vlc_demodulate(demod_data, demodulator, qam_handle);
    
    %% Decoding
    if FEC_Coding ~= 0
        
        deleaving = @de_interleaving;
        deleav_data = vlc_decode(demod_data, deleaving, interleave_table);
        
        decoder = @conv_decoder;
        decoder_data= vlc_decode(deleav_data, decoder, CodeTrellis);
    else
        decoder_data = demod_data;
    end
    
    %% System analysis
    [number_of_errors,bit_error_rate] = biterr(in_data,decoder_data);
    BER(k) = bit_error_rate;
end

%BER

%% Plot
figure;
semilogy(SNR_elec,BER,'r');

figure;
semilogy(SNR_elec_nodc,BER,'g');

figure;
semilogy(SNR_opt,BER,'b');
