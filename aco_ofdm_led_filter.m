function out = aco_ofdm_led_filter(in, led)
%clipping the negative part
%clipping 
out = in;
out(out<led.min) = led.min;
out(out>led.max) = led.max;