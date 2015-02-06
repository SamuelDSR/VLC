function out = qam_demodulator(in, qamDemodParams)
    bitAlloc    = cell2mat(qamDemodParams(1));
    powerAlloc  = cell2mat(qamDemodParams(2));
    
    nSubcar             = length(bitAlloc);   % number of nSubcar per ofdm symbol
    nBitsPerOfdmSymbol  = sum(bitAlloc);      % number of bits per ofdm symbol
    nOfdmSymbol         = length(in)/nSubcar; % number of ofdm symbol
    
    in              = reshape(in,nSubcar,nOfdmSymbol); 
    out             = zeros(nBitsPerOfdmSymbol,nOfdmSymbol);
    
    bitsIndex       = 0;
    
    % de bit loading, different nSubcar has different constellation size
    for i = 1:nSubcar
        
        if bitAlloc(i) ~= 0 % if the subcarrier is modulated, demodulated with qam, else do nothing
            
            order = 2^bitAlloc(i);
            
            handle = modem.qamdemod('M',order,'SymbolOrder','gray','OutputType','bit','DecisionType','hard decision');
            
            normalizefactor = 2/3*(order-1);
            
            out(bitsIndex+1:bitsIndex+bitAlloc(i),:) = demodulate(handle,sqrt(normalizefactor)/powerAlloc(i)*in(i,:));
            
            %             % de power loading
            %             if powerAlloc(i) ~= 0 %if the subcarrier has power allocated
            %                 out(bitsIndex+1:bitsIndex+bitAlloc(i),:) = demodulate(handle,sqrt(normalizefactor)/powerAlloc(i)*in(i,:));
            %             else
            %                 out(bitsIndex+1:bitsIndex+bitAlloc(i),:) = 0;
            %             end
            
            bitsIndex = bitsIndex+bitAlloc(i);     
            
        end

    end
    
    out  = reshape(out,nBitsPerOfdmSymbol*nOfdmSymbol,1);
    
end