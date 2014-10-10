function out = qam_demodulator(in, handle)
    deModHandle = modem.qamdemod(handle);
    out = demodulate(deModHandle,in);
end