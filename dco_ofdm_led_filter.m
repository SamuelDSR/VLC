function out = dco_ofdm_led_filter(in, led)
%get dc bias
out = in + led.dc_bias;
%clipping 
out(out<led.min) = led.min;
out(out>led.max) = led.max;
