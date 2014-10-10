%% generate random signal
Nsym = 10000;
nBitPerSymbol = 4;
FEC_Coding = 0;

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
ModOrder = 16;
modulator = @qam_modulator;
[mod_data, qam_handle] = vlc_modulate(enc_data, modulator, ModOrder);

%CAP modulator
FilterSpan = 8;
UpSamplingFactor = 4;
modulator_param = [FilterSpan,UpSamplingFactor];
modulator = @cap_modulator;
[mod_data, cap_filters] = vlc_modulate(mod_data, modulator, modulator_param);

%% LED filtering
led = [];
tx_data  = vlc_led_filter(mod_data, led);

%% channel filtering
channel = [];
rx_data  = vlc_channel_filter(tx_data, channel);

%% detection using receiver filter
pd = [];
detect_data = vlc_pd_filter(rx_data, pd);

%% Add gaussian noise
EbNo = 30;
SNR = EbNo + 10*log10(nBitPerSymbol) - 10*log10(UpSamplingFactor);
detect_data = awgn(detect_data,SNR,'measured');

%% Demodulation
demodulator = @cap_demodulator;
demod_data  = vlc_demodulate(detect_data, demodulator, cap_filters);

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
[number_of_errors,bit_error_rate] = biterr(source,decoder_data)
