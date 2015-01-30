% generate Ac+dc signal
f = 100; %100Hz, frequency of ac component
fs = 250; %sampling frequency
t = 0:1/fs:1.5;
ac = sin(2*pi*t);
dc = 0.315;
signal = ac + dc;

plot(t,signal);

% creat a high pass digital filter
