function out = qam_demodulator(in, qamDemodParams)
    bitAlloc    = cell2mat(qamDemodParams(1));
    powerAlloc  = cell2mat(qamDemodParams(2));
    
    nSubcar         = length(bitAlloc);   % number of nSubcar per ofdm symbol
    nBitsPerSymbol  = sum(bitAlloc);      % number of bits per ofdm symbol
    nSymbol         = length(in)/nSubcar; % number of ofdm symbol
    
    in              = reshape(in,nSubcar,nSymbol); 
    out             = zeros(nBitsPerSymbol,nSymbol);
    
    bitsIndex       = 0;
    
    % de bit loading, different nSubcar has different constellation size
    for i = 1:nSubcar
        
        order = 2^bitAlloc(i);
        
        handle = modem.qamdemod('M',order,'SymbolOrder','gray','OutputType','bit','DecisionType','hard decision');
            
        normalizefactor = 2/3*(order-1);
        
        % de power loading
        out(bitsIndex+1:bitsIndex+bitAlloc(i),:) = demodulate(handle,sqrt(normalizefactor)/powerAlloc(i)*in(i,:));
        
        bitsIndex = bitsIndex+bitAlloc(i);
    end
    
    out  = reshape(out,nBitsPerSymbol*nSymbol,1);
end