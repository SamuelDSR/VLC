function [out, blkSize] = dco_ofdm_modulator(in, modulator_params)
% get parameters
nSubcar     = modulator_params(1);
cpSize      = modulator_params(2);
fftSize     = (nSubcar+1)*2;

if rem(length(in),nSubcar) == 0
    
    nOfdmSymbol = length(in)/nSubcar;   
    blkSize     = fftSize+cpSize;
    out         = zeros(nOfdmSymbol*blkSize,1);  %pre-allocate memeory
    
    for i = 1:nOfdmSymbol
  
        in_temps        = in(1+(i-1)*nSubcar:i*nSubcar);
        
        % pilot insertion and hermetian sysmetry to ensure a real data after FFT
        pilot_ins_data  = [0; in_temps;  0; (fliplr(conj(in_temps)'))'];

        % fourier transform time domain data (note that we multiply a factor of sqrt(FFT Size) to make sure we have same average energy after IFFF)
        IFFT_data       = sqrt(fftSize)*ifft(pilot_ins_data);
        
        % add cycle prefix 
        out((i-1)*blkSize+1:i*blkSize) = [IFFT_data(end-cpSize+1:end); IFFT_data];
    end
else
    error('DCO-OFDM wrong input size for modulation');
end
