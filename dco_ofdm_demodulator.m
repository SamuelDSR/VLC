function out = dco_ofdm_demodulator(in, demod_params)
%get params
subcar   = demod_params(1);
cp_size  = demod_params(2);
fft_size = (subcar+1)*2; % FFT size (this relation is only true in ACO-OFDM)
blk_size = fft_size+cp_size;  % OFDM symbol length

if mod(length(in),blk_size) == 0  
    
    Nsym = length(in)/blk_size;
    
    out = zeros(Nsym*subcar,1);  %pre-allocate memeory
    
    for i = 1:Nsym
  
        in_blk = in(1+(i-1)*blk_size:i*blk_size);
        
        % cp removal
        in_rmcp = in_blk(cp_size+1:end);
        
        % fft
        in_fft = fft(in_rmcp);
        
        % pilot removal for ACO-OFDM
        in_subcar = in_fft(2:subcar+1);
        
        out(1+(i-1)*subcar:i*subcar) = in_subcar;       
    end
    
else
    error('ACO-OFDM wrong input size for demodulation');
end