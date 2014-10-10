function out = conv_decoder(in, t)
    tblen = log2(t.numOutputSymbols);
    out = vitdec(in, t, tblen, 'cont', 'hard');
end