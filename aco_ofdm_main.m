close all;
clear all;

%% Parameters
%using FEC or not
FEC_Coding = 0;

%SNR 
N0 = 10^-21;  % spectral density level of gaussian noise at the receiver
B = 20*10^6;  % signal bandwidth B = 20MHz

SNR = 100:5:140; %dB, Desired electrical SNR per bit before clipping at the transmitter (excluding the dc_bias)


%ACO-OFDM params
Nsym = 10000;  %number of OFDM symbols
subcar = 32;
FFT_size = subcar*4;
cp_size = 16;

nBitPerSymbol = 4;
yita = nBitPerSymbol*subcar/(4*subcar+cp_size); % Spectral efficiency

%LED params
led.min = 0.1; % A
led.max = 1 ; % A
led.dc_bias = led.min; %set the dc bias according to the signal variance

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
PAPR_dB = zeros(size(SNR));

for k = 1:length(SNR)
    
    %% generate random signal
    generator = @random_bin_generator;
    in_data =  vlc_source_generate(Nsym*nBitPerSymbol*subcar, generator);
    
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
    
    %ACO-OFDM modulator
    aco_params = [subcar,cp_size];
    modulator = @aco_ofdm_modulator;
    [mod_data, blk_size] = vlc_modulate(mod_data, modulator, aco_params);
    
    %% Scaling the signal for matching the desired SNR
    P_elec = B*yita*N0*10^(SNR(k)/10);
    signal_amp_moyen = sum(mod_data)/length(mod_data);
    signal_var = sum(mod_data.^2)/length(mod_data);
    signal_peak = max(mod_data.^2);
    PAPR_dB(k) = 10*log10(signal_peak/signal_var);
    
    alpha = sqrt(P_elec/signal_var); %scaling factor for fitting to the LED modulation interval
    mod_data = alpha*mod_data;
    
    %% LED Clipping  
    led_clipping=@aco_ofdm_led_filter;
    tx_data  = vlc_led_filter(mod_data, led_clipping, led);
    
    %% Calculating the Average P_elec and Average P_opt after LED filter
    P_elec_avg = sum(tx_data(1:1000*blk_size).^2)/blk_size/1000;  %signal elec after clipping including dc-bias
    P_elec_avg_nodc = sum((tx_data(1:1000*blk_size)-led.dc_bias).^2)/blk_size/1000; %signal elec after clipping excluding dc-bias
    
    P_opt_avg = sum(tx_data(1:1000*blk_size))/blk_size/1000;
    
    SNR_elec(k) = 10*log10(P_elec_avg/(N0*B*yita));
    SNR_elec_nodc(k) = 10*log10(P_elec_avg_nodc/(N0*B*yita));
    SNR_opt(k) = 10*log10(P_opt_avg/(N0*B*yita));
    
    %figure;
    %plot(tx_data(1:100*blk_size));
    
    %% channel filtering
    rx_data  = vlc_channel_filter(tx_data, channel);
    
    %% detection using receiver filter
    detect_data = vlc_pd_filter(rx_data, pd);
    
    %% Add gaussian noise
    noise_var = B*N0;
    detect_data = detect_data + sqrt(noise_var)*randn(size(detect_data));
    
    %% Demodulation
    pilot_info = alpha*channel*pd;
    aco_params(3) = pilot_info;
    
    demodulator = @aco_ofdm_demodulator;
    demod_data  = vlc_demodulate(detect_data, demodulator, aco_params);
    
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

% BER
% 
% SNR_elec
% SNR_elec_nodc
% SNR_opt
%% Plot
figure;
semilogy(SNR_elec,BER,'r');

figure;
semilogy(SNR_elec_nodc,BER,'g');

figure;
semilogy(SNR_opt,BER,'b');

