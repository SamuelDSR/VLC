function [out, blk_size] = dco_ofdm_modulator(in, modulator_params)
% get parameters
subcar   = modulator_params(1);
cp_size  = modulator_params(2);
fft_size = (subcar+1)*2;

if mod(length(in),subcar) == 0  
    
    Nsym = length(in)/subcar;
    
    blk_size = fft_size+cp_size;
    
    out = zeros(Nsym*blk_size,1);  %pre-allocate memeory
    
    for i = 1:Nsym
  
        in_temps = in(1+(i-1)*subcar:i*subcar);
        
        % pilot insertion and hermetian sysmetry to ensure a real data after FFT
        pilot_ins_data = [0; in_temps;  0; (fliplr(conj(in_temps)'))'];

        % fourier transform time doamain data
        IFFT_data =ifft(pilot_ins_data);
        
        % add cycle prefix 
        out((i-1)*blk_size+1:i*blk_size) = [IFFT_data(end-cp_size+1:end); IFFT_data];
    end
else
    error('DCO-OFDM wrong input size for modulation');
end
