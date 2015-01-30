clear all;
close all;
ledCutoffFrequecy = 5*10^6;
ledFilterOrder    = 31;
B                 = 100*10^6; 
ledFilterCoeff    = led_lp_channel(ledFilterOrder,ledCutoffFrequecy/B,1);

data = rand(30,1);
txData = [zeros(10,1);data];
txData(1:10) = txData(end-10+1:end);

%txData = repmat(tx_Data,20,1);
rxData = vlc_channel_filter(txData,ledFilterCoeff);
rxData = rxData(1:30);

plot(data,'r');
hold on
plot(rxData,'b');

Hf_channel = abs(fft(ledFilterCoeff,32));
Hf_txData  = abs(fft(data,32));
Hf_rxData  = abs(fft(rxData,32));

figure
plot(Hf_txData,'r');
hold on
plot(Hf_rxData./Hf_channel,'b');

err = data - ifft(Hf_rxData./Hf_channel,30)