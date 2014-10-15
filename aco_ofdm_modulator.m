function [out, blk_size] = aco_ofdm_modulator(in, modulator_param)
% get parameters
subcar   = modulator_param(1);
cp_size  = modulator_param(2);

fft_size = (subcar+1)*4; % it's in ACO-OFDM

if mod(length(in),subcar) == 0  %ensure the input data length is a multiple of fft_size
    
    Nsym = length(in)/subcar;
    
    blk_size = fft_size+cp_size;
    
    out = zeros(Nsym*blk_size,1);  %pre-allocate memeory
    
    for i = 1:Nsym
  
        in_temps = in(1+(i-1)*subcar:i*subcar);
        
        % pilot insertion and hermetian sysmetry to ensure a real data after FFT
        pilot_ins_data = zeros(fft_size,1);
        for j = 1:subcar
            pilot_ins_data(j*2) = in_temps(j);
            pilot_ins_data(fft_size - j*2 + 2) = conj(in_temps(j));
        end

        % fourier transform time doamain data
        IFFT_data =ifft(pilot_ins_data);
        
        %clipping
        IFFT_data(IFFT_data<0)=0; % clipping the negative part

        % add cycle prefix 
        out((i-1)*blk_size+1:i*blk_size) = [IFFT_data(end-cp_size+1:end); IFFT_data];
    end
else
    error('ACO-OFDM wrong input size for modulation');
end
        
        
        