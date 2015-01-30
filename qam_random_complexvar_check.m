clear all;
close all;
%% power bit loading checking
subcar = 127;
bitspersymbol = randi([2 10],subcar,1);

powerloading = rand(128,1);
%powerloading = 1/subcar*ones(subcar,1);
powerloading = powerloading/sum(powerloading);

out = zeros(2000,subcar);
symbols = 1000;

for i=1:subcar    
    handle = modem.qammod(2^bitspersymbol(i));
    handle.InputType = 'Bit';
    handle.SymbolOrder = 'Gray';
    in = randi([0 1],bitspersymbol(i)*2000,1);
    normalizefactor = 2/3*(2^bitspersymbol(i)-1);
    out(:,i) = modulate(handle,in);
    out(:,i) = 1/normalizefactor*powerloading(i)*out(:,1);
end

%% FFT (DC0-OFDM)
out = [zeros(2000,1),out,zeros(2000,1),(fliplr(conj(out)))];
fft_res = 1/sqrt(subcar*2+2)*fft(out',subcar*2+2)';

%% Gaussian variable with mean zero checking
mean = zeros(2000,1);
for i=1:2000
    for j = 1:subcar*2+2
        mean(i) = mean(i) + fft_res(i,j);
    end
    mean(i) = mean(i)/(subcar*2+2);
end

%% Normality test
test = zeros(2000,1);
for i = 1:2000
    case1 = real(fft_res(i,:));
    var = sum(case1.^2)/length(case1);
    case1 = case1/sqrt(var);
    test(i) = kstest(case1);
end

    