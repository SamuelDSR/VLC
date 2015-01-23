subcar = 16;
order = [10 10 9 8 8 7 7 6 6 6 5 4 3 3 2 1];
out = zeros(2000,subcar);

for i=1:subcar    
    handle = modem.qammod(2^order(i));
    handle.InputType = 'Bit';
    handle.SymbolOrder = 'Gray';
    in = randi([0 1],order(i)*2000,1);
    out(:,i) = modulate(handle,in);
end