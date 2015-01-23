function out = aco_ofdm_demodulator(in, demod_params)

subcar   = demod_params(1);
cp_size  = demod_params(2);
attenuation = demod_params(3);

fft_size = subcar*4; % FFT size (this relation is only true in ACO-OFDM)
blk_size = fft_size+cp_size;  % OFDM symbol length

if mod(length(in),blk_size) == 0  
    
    Nsym = length(in)/blk_size;
    
    blk_size = fft_size+cp_size;
    
    out = zeros(Nsym*subcar,1);  %pre-allocate memeory
    
    for i = 1:Nsym
  
        in_blk = in(1+(i-1)*blk_size:i*blk_size);
        
        % cp removal
        in_rmcp = in_blk(cp_size+1:end);
        
        % composenate the attenuation
        in_rmcp = in_rmcp*2/attenuation;
        
        % fft
        in_fft = fft(in_rmcp);
        
        % pilot removal for ACO-OFDM
        in_subcar = zeros(subcar,1);
        
        for j = 1:subcar
            in_subcar(j,1) = in_fft(j*2);
        end
        
        out(1+(i-1)*subcar:i*subcar) = in_subcar;       
    end
    
else
    error('ACO-OFDM wrong input size for demodulation');
end