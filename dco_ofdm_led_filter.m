function out = dco_ofdm_led_filter(in, led)
%get dc bias
out = in + led.dcBias;
%clipping 
out(out<led.minCurrent) = led.minCurrent;
out(out>led.maxCurrent) = led.maxCurrent;
