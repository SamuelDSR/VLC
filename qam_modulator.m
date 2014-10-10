function [out, handle] = qam_modulator(in, order)
    handle = modem.qammod(order);
    handle.InputType = 'Bit';
    handle.SymbolOrder = 'Gray';
    out = modulate(handle,in);
end