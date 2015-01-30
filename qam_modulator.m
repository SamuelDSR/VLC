% Qam modulation using generated bit and power loading schemes
% =============================================================
% Input Params:
%
% in: 1D-array input binary data

% qamParams: should contain 2 params:
%   1. bitAlloc:  bit loading schemes
%   2. energyAlloc: power loading schemes
%
% =============================================================
% Output params:
%
% out: 1D-array of length nSubcar*nOfdmSymbol
% remBits:  rem bits not modulated
function [out, remBits] = qam_modulator(in, qamParams)

bitAlloc = cell2mat(qamParams(1));
powerAlloc = cell2mat(qamParams(2));
nSubcar = length(bitAlloc);

totalBits = length(in);                     % total bits need to modulate
nBitsPerSymbol = sum(bitAlloc);             % number of bits per ofdm symbol
nOfdmSymbol = floor(totalBits/nBitsPerSymbol);  % number of ofdm symbol, floor rounding
remBits = rem(totalBits, nBitsPerSymbol);   % rem bits not modulated (length of rem bits not long enough to fill in a ofdm symbol)


in = reshape(in(1:nOfdmSymbol*nBitsPerSymbol),nBitsPerSymbol,nOfdmSymbol);
out= zeros(nSubcar,nOfdmSymbol);

bitsIndex = 0;
for i = 1:length(bitAlloc)
    order = 2^bitAlloc(i);
    handle = modem.qammod(order);
    handle.InputType = 'Bit';
    handle.SymbolOrder = 'Gray';     
    normalizefactor = 2/3*(order-1);
    out(i,:) = sqrt(powerAlloc(i))/sqrt(normalizefactor)*modulate(handle,in(bitsIndex+1:bitsIndex+bitAlloc(i),:));
    bitsIndex = bitsIndex+bitAlloc(i);
end

out = reshape(out,nSubcar*nOfdmSymbol,1);  %% note that the function of  reshape is taken columwise, so the symbol order are not corrupted

