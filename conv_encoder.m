function [out, CodeTrellis] = conv_encoder(in, encoder_param)
    CodeTrellis = poly2trellis([5 4],[23 35 0;0 5 13]);
    if mod(length(in), CodeTrellis.numOutputSymbols) == 0
        out = convenc(in, CodeTrellis);
    else
        error('Wrong input length to convolution encoder');
    end
        