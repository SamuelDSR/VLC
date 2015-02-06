  function [shanonCapacity powerAllocated] = ofdm_waterfilling(nSubChannel,totalPower,channelStateInformation,bandwidth,noiseDensity)
%==========================================================================
% by:
%   Hamid Ramezani
%   my email : hamid.ramezani@gamil.com
%
% Date:
%   20 OCT 2006 , version 1
%
%
% to maximize the capacity of a frequency selective channel, the
% waterfilling algorithm is used. OFDM modulation devides the total
% bandwidth into N subchannels. Large number of subchannels cause each
% subchannel experience a flat fading channel if the proper cyclic prefix
% is added to the end of each OFDM symbol. The waterfilling algorithm
% assigns more power to subchannels which experience good condition and may
% assign no power to bad conditioned subchannels (subchannels withdeep fading).
%
%   for further information see:
%     J. A. C. Bingham, “Multicarrier Modulation for Data Transmission:
%     an Idea whose Time Has Come,” IEEE Communications Magazine, vol. 28,
%     no. 5, pp. 5-14, May 1990
%
% =========================================================================
%                      Parameter definition
% =========================================================================
%
% nSubChannel               : Number of subchannels (4,16,32,...,2^N)
% totalPower                : Total available power for each OFDM symbol
%                             (p watt)
% channelStateInformation   : Channel state information 1_by_nSubChannel
%                             matrix "random('rayleigh',1,1,nSubChannel)"
% bandwidth                 : Total available bandwidth (in Hz)
% noiseDensity              : one sided noise spectral dencity (watt/Hz)

% =========================================================================
%                               Example
% =========================================================================
% nSubChannel             = 16;
% totalPower              = 1e-5;           %  -20 dBm
% channelStateInformation = random('rayleigh',1/0.6552,1,nSubChannel);
% bandwidth               = 1e6;            %  1 MHz
% noiseDensity            = 1e-11;          % -80 dBm
% [Capacity PowerAllocated] = ofdmwaterfilling(...
%    nSubChannel,totalPower,channelStateInformation,bandwidth,noiseDensity)



%< Parameter Computation >

subchannelNoise     = ...
    noiseDensity*bandwidth/nSubChannel;
carrierToNoiseRatio = ...
    channelStateInformation.^2/subchannelNoise;

initPowerAllo       = ...                   (Formula 1)
    (totalPower + sum(1./carrierToNoiseRatio))...
    /nSubChannel - 1./carrierToNoiseRatio;

% Until this stage maybe some subchannels were allocated negetive power
% becouse of the optimization formula 1. At the following stage these
% subchannels will be eliminated from using any power and the algorithm is
% repeated for the other subchannels, until all the remained ones, will
% be allocated the positive power.

% < Iterative part of the  algorithm >
while(length( find(initPowerAllo < 0 )) > 0 )
    negIndex       = find(initPowerAllo <= 0);
    posIndex       = find(initPowerAllo >  0);
    nSubchannelRem = length(posIndex);
    initPowerAllo(negIndex) = 0;
    CnrRem         = carrierToNoiseRatio(posIndex);
    powerAlloTemp  = (totalPower + sum(1./CnrRem))...  
        /nSubchannelRem - 1./CnrRem;
    initPowerAllo(posIndex) = powerAlloTemp;
end

% < Output Computation >

% amount of power allocated to each subchannel
powerAllocated = initPowerAllo';
% total capacity of a channel based on shanon theory
shanonCapacity = bandwidth/nSubChannel * ...
    sum(log2(1 + initPowerAllo.*carrierToNoiseRatio));

% <Graphical Observation>

%By observing the figure, it is clear that the power like water fills the
%container which is made by noise to carrier ratio or channel state
%information

f1 = figure(2);
clf;
set(f1,'Color',[1 1 1]);
bar((initPowerAllo + 1./carrierToNoiseRatio),1,'r')
hold on;
bar(1./carrierToNoiseRatio,1);
xlabel('subchannel indices');
title('Water filling algorithm')

legend('amount of power allocated to each subchannel',...
    'Noise to Carrier Ratio')

