function out = cap_demodulator(in, cap_filters)

    inphase_pulse = cell2mat(cap_filters(1));
    quadrature_pulse = cell2mat(cap_filters(2));
    FilterTaps =  length(inphase_pulse)-1;
    UpSamplingFactor = cell2mat(cap_filters(3));
    
    % calculated filter gain
    inphase_inphase_pulse = conv(fliplr(inphase_pulse),inphase_pulse);
    quadrature_quadrature_pulse = conv(fliplr(quadrature_pulse),quadrature_pulse);
    %inphase_quadrature_pulse=conv(fliplr(quadrature_pulse),inphase_pulse);
    inphase_gain = inphase_inphase_pulse(FilterTaps+1);
    quadrature_gain = quadrature_quadrature_pulse(FilterTaps+1);
    
    %Filter the received signal with the corresponding MF filter
    in(end+1:end+FilterTaps/2) = 0; %zero padding

    inphaseRx = filter(fliplr(inphase_pulse),1,in);
    quadratureRx = filter(fliplr(quadrature_pulse),1,in);

    %Delay
    inphaseRx = inphaseRx(FilterTaps/2+1:end);
    quadratureRx = quadratureRx(FilterTaps/2+1:end);
    
    %Downsampling 
    inphaseRx = downsample(inphaseRx,UpSamplingFactor);
    quadratureRx = downsample(quadratureRx, UpSamplingFactor);

    %Normalization
    inphaseRx = inphaseRx/inphase_gain;
    quadratureRx = quadratureRx/quadrature_gain;

    % Summation
    out = inphaseRx - quadratureRx*1i;