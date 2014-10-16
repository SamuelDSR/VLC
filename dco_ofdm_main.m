%% generate random signal
%number of OFDM symbols
Nsym = 10000;  
%aco ofdm params
subcar = 15;
cp_size = 3;
nBitPerSymbol = 4;

%using FEC or not
FEC_Coding = 0;

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
dco_params = [subcar,cp_size];
modulator = @dco_ofdm_modulator;
[mod_data, blk_size] = vlc_modulate(mod_data, modulator, dco_params);

%% LED filtering
led.dc_bias = 2;
led.min = 0;
led.max = 100;
led_filter=@dco_ofdm_led_filter;
tx_data  = vlc_led_filter(mod_data, led_filter, led);

%% channel filtering
channel = [];
rx_data  = vlc_channel_filter(tx_data, channel);

%% detection using receiver filter
pd = [];
detect_data = vlc_pd_filter(rx_data, pd);

%% Add gaussian noise

EbNo = 20; 
snr = EbNo + 10*log10(subcar/blk_size) + 10*log10(log2(ModOrder));
detect_data = awgn(detect_data,snr,'measured');

%% Demodulation
demodulator = @dco_ofdm_demodulator;
demod_data  = vlc_demodulate(detect_data, demodulator, dco_params);

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
[number_of_errors,bit_error_rate] = biterr(in_data,decoder_data)
