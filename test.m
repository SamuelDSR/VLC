clear all;
close all;

ledCutoffFrequecy = 5*10^6;
ledFilterOrder    = 100;
B                 = 200*10^6; 
ledFilterCoeff    = led_lp_channel(ledFilterOrder,ledCutoffFrequecy/B,1);

fvtool(ledFilterCoeff,1)
%    = ledFilterCoeff + 1i*ledFilterCoeff;
% 
% ledFilterCoeff2   = led_lp_channel(ledFilterOrder*2,ledCutoffFrequecy/B/2,1);
% plot(ledFilterCoeff1);
% %hold on
% figure
% plot(1:length(ledFilterCoeff1),ledFilterCoeff2(1:2:length(ledFilterCoeff2)))

 data_f = 0.001*randi([2,5],256,1);
 data = ifft(data_f);

cp = 36;
%data = randi([2,4],112,1);
txData = [zeros(cp,1);data];
txData(1:cp) = data(256-cp+1:256);

%txData = repmat(txData,20,1);
rxData = vlc_channel_filter(txData,ledFilterCoeff);

% plot(data,'r');   
% hold on
% plot(rxData,'b');
rxData = rxData(cp+1:256+cp);
rxData_f = fft(rxData,256);

err = rxData_f./fft(ledFilterCoeff,256) - data_f;
plot(10*log(abs(err)));
grid on
% Hf_channel = abs(fft(ledFilterCoeff,128));
% Hf_txData  = abs(fft(data,128));
% Hf_rxData  = abs(fft(rxData(16+1:112+16),128));
% 
% figure
% plot(Hf_txData,'r');
% hold on
% plot(Hf_rxData./Hf_channel,'b');
% 
% 
% %err = data - ifft(Hf_rxData./Hf_channel,30)