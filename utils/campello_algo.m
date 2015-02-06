% campello_algo.m
% This function is used by Campello's algorithm to allocate bits and energy for 
% each subchannel optimally.
% Input:
%       bitAlloc    : initial bit allocation scheme
%       gap         : channel gap
%       gainToNoise : channel gain of sub channel versus Noise 
%       nSubcar     : total Number of sub channel
%       maxOrder    : max constellation size

function [bitAlloc, powerAlloc] = campello_algo(bitAlloc,gap,gainToNoise,nSubcar,totalPower,maxOrder)

%% calculation of power table
% add two rows (bit=0,Energy=0) and(bit=maxOrder+1,Energy=Inf) for convience in programming
powerTable = zeros(maxOrder+2,nSubcar);
for i = 1:nSubcar
    for j = 2:maxOrder+1
        powerTable(j,i) = 2*gap*(2^j-1)/gainToNoise(i);
    end
end
powerTable(maxOrder+2,:) = Inf;

%% calculation of incremental energy table

increPowerTable = zeros(maxOrder+2,nSubcar);  
for i = 2:maxOrder+2
    increPowerTable(i,:) = powerTable(i,:) - powerTable(i-1,:);
end


%% current array of energy gained to decrease one bit
bitDecreaseTable = zeros(1,nSubcar);
for i = 1:nSubcar
    bitDecreaseTable(i) = increPowerTable(bitAlloc(i)+1,i);  %here, plus 1 because in matlab array indice begin with 1
end
% max energy gained to decrease one bit and the correspondant subcar indice
[maxDecreasePower, maxIndice] = max(bitDecreaseTable);

%% current array of energy needed to increase one bit
bitIncreaseTable = zeros(1,nSubcar);
for i = 1:nSubcar
    bitIncreaseTable(i) = increPowerTable(bitAlloc(i)+2,i);
end
% min energy needed to increase one bit and the correspondant subcar indice
[minIncreasePower, minIndice] = min(bitIncreaseTable);

%% step of efficientizing
while minIncreasePower < maxDecreasePower
    % change bit loading scheme
    bitAlloc(minIndice) = bitAlloc(minIndice)+1;
    bitAlloc(maxIndice) = bitAlloc(maxIndice)-1;
    % updating
    for i = 1:nSubcar
        bitDecreaseTable(i) = increPowerTable(bitAlloc(i)+1,i);
    end
    [maxDecreasePower, maxIndice] = max(bitDecreaseTable);
    
    bitIncreaseTable = zeros(1,nSubcar);
    for i = 1:nSubcar
        bitIncreaseTable(i) = increPowerTable(bitAlloc(i)+2,i);
    end
    [minIncreasePower, minIndice] = min(bitIncreaseTable);
    %     minIncreasePower
    %     maxDecreasePower
end

%% current Power allocation table
currentPower  = 0;
for i=1:nSubcar
    currentPower = currentPower+powerTable(bitAlloc(i)+1,i);
end

%% step of tighting (to satisfait the total energy constraint)

while  totalPower-currentPower < 0 || totalPower-currentPower >= minIncreasePower
    
    if totalPower-currentPower < 0 % energy exceeding limit
        currentPower = currentPower - maxDecreasePower;
        bitAlloc(maxIndice) = bitAlloc(maxIndice)-1;
        % updating
        for i = 1:nSubcar
            bitDecreaseTable(i) = increPowerTable(bitAlloc(i)+1,i);
        end
        [maxDecreasePower, maxIndice] = max(bitDecreaseTable);
        
        bitIncreaseTable = zeros(1,nSubcar);
        for i = 1:nSubcar
            bitIncreaseTable(i) = increPowerTable(bitAlloc(i)+2,i);
        end
        [minIncreasePower, minIndice] = min(bitIncreaseTable);
        
    else
        currentPower = currentPower + minIncreasePower;
        bitAlloc(minIndice) = bitAlloc(minIndice)+1;
        % updating
        for i = 1:nSubcar
            bitDecreaseTable(i) = increPowerTable(bitAlloc(i)+1,i);
        end
        [maxDecreasePower, maxIndice] = max(bitDecreaseTable);
        
        bitIncreaseTable = zeros(1,nSubcar);
        for i = 1:nSubcar
            bitIncreaseTable(i) = increPowerTable(bitAlloc(i)+2,i);
        end
        [minIncreasePower, minIndice] = min(bitIncreaseTable);
    end
end

% get the energy allocation scheme
powerAlloc = zeros(nSubcar,1);
for  i=1:nSubcar
    powerAlloc(i) = powerTable(bitAlloc(i)+1,i);
end

%% assign the left energy

leftPower1 = totalPower - currentPower;
leftPower2 = totalPower - sum(powerAlloc);
leftPower2

if abs(leftPower1-leftPower2)<1e-10  % assign the additional power evenly to all modulated sub carrier
    nModSubcar = find(bitAlloc~=0);
    powerAlloc(nModSubcar) = powerAlloc(nModSubcar) + leftPower1/length(nModSubcar);
else
    error('Something is wrong about the bit-power loading algorithme, check the program');
end
    
%% visualization of power-bit loading
figure
bar(1:nSubcar,powerAlloc+1./gainToNoise,1,'r')
hold on;
bar(1./gainToNoise,1);
xlabel('subchannel indices');
title('Power loading')
legend('amount of power allocated to each subchannel',...
    'Noise to Carrier Gain Ratio')

figure
bar(1:nSubcar,bitAlloc);
xlabel('subchannel indices');
title('Bit loading')
legend('constellation size of each subchannel')


