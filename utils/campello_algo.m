% campello_algo.m
% This function is used by Campello's algorithm to allocate bits and energy for 
% each subchannel optimally.
% Input:
%       bitAlloc    : initial bit allocation scheme
%       powerAlloc  : initital power allocation scheme
%       gap         : channel gap
%       gainToNoise : channel gain of sub channel versus Noise 
%       nSubcar     : total Number of sub channel
%       maxOrder    : max constellation size

function [bitAlloc, energyAlloc] = campello_algo(bitAlloc,energyAlloc,gap,gainToNoise,nSubcar,totalEnergy,maxOrder)

%% calculation of energy table
% add two rows (bit=0,Energy=0) and(bit=maxOrder+1,Energy=Inf) for convience in programming
energyTable = zeros(maxOrder+2,nSubcar);
for i = 1:nSubcar
    for j = 2:maxOrder+1
        energyTable(i,j) = 2*gap*gainToNoise(i)*(2^j-1);
    end
end
energyTable(maxOrder+2:i) = Inf;

%% calculationb of incremental energy table

increEnergyTable = zeros(maxOrder+2,nSubcar);  
for i = 2:maxOrder+2
    increEnergyTable(i,:) = energyTable(i,:) - energyTable(i-1,:);
end


%% current array of energy gained to decrease one bit
bitDecreaseTable = zeros(1,nSubcar);
for i = 1:nSubcar
    bitDecreaseTable(i) = increEnergyTable(bitAlloc(i)+1,i);  %here, plus 1 because in matlab array indice begin with 1
end
% max energy gained to decrease one bit and the correspondant subcar indice
[maxDecreaseEnergy, maxIndice] = min(bitDecreaseTable);

%% current array of energy needed to increase one bit
bitIncreaseTable = zeros(1,nSubcar);
for i = 1:nSubcar
    bitIncreaseTable(i) = increEnergyTable(bitAlloc(i)+1,i);
end
% min energy needed to increase one bit and the correspondant subcar indice
[minIncreaseEnergy, minIndice] = min(bitIncreaseTable);

%% step of efficientizing
while minIncreaseEnergy < maxDecreaseEnergy
    % change bit loading scheme
    bitAlloc(minIndice) = bitAlloc(minIndice)+1;
    bitAlloc(maxIndice) = bitAlloc(maxIndice)-1;
    % updating
    for i = 1:nSubcar
        bitDecreaseTable(i) = increEnergyTable(bitAlloc(i)+1,i);
    end
    [maxDecreaseEnergy, maxIndice] = min(bitDecreaseTable);
    
    bitIncreaseTable = zeros(1,nSubcar);
    for i = 1:nSubcar
        bitIncreaseTable(i) = increEnergyTable(bitAlloc(i)+1,i);
    end
    [minIncreaseEnergy, minIndice] = min(bitIncreaseTable);
end

%% current Energy allocation table
currentEnergy  = 0;
for i=1:nSubcar
    currentEnergy = currentEnergy+energyTable(bitAlloc(i),i);
end

%% step of tighting (to satisfait the total energy constraint)

while  currentEnergy < totalEnergy || currentEnergy >= minIncreaseEnergy
    
    if currentEnergy < totalEnergy % energy exceeding limit
        currentEnergy = currentEnergy - maxDecreaseEnergy;
        bitAlloc(maxIndice) = bitAlloc(maxIndice)-1;
        % updating
        for i = 1:nSubcar
            bitDecreaseTable(i) = increEnergyTable(bitAlloc(i)+1,i);
        end
        [maxDecreaseEnergy, maxIndice] = min(bitDecreaseTable);
        
        bitIncreaseTable = zeros(1,nSubcar);
        for i = 1:nSubcar
            bitIncreaseTable(i) = increEnergyTable(bitAlloc(i)+1,i);
        end
        [minIncreaseEnergy, minIndice] = min(bitIncreaseTable);
        
    else
        currentEnergy = currentEnergy + minIncreaseEnergy;
        bitAlloc(minIndice) = bitAlloc(minIndice)+1;
        % updating
        for i = 1:nSubcar
            bitDecreaseTable(i) = increEnergyTable(bitAlloc(i)+1,i);
        end
        [maxDecreaseEnergy, maxIndice] = min(bitDecreaseTable);
        
        bitIncreaseTable = zeros(1,nSubcar);
        for i = 1:nSubcar
            bitIncreaseTable(i) = increEnergyTable(bitAlloc(i)+1,i);
        end
        [minIncreaseEnergy, minIndice] = min(bitIncreaseTable);
    end
end

% get the energy allocation scheme

for  i=1:nSubcar
    energyAlloc = energyTable(bitAlloc(i)+1,i);
end

%% assign the left energy

leftEnergy1 = totalEnergy - currentEnergy;
leftEnergy2 = totalEnergy - sum(energyAlloc);
leftEnergy2

if leftEnergy1 == leftEnergy  % assign the additional power evenly to all subcar
    energyAlloc = energyAlloc + leftEnergy1/nSubcar;  
else
    error('Something is wrong about the bit-power loading algorithme, check the program');
end
    
