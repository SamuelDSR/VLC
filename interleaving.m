function [int_aft_interleav, properties] = interleaving( int_pre_interleav ,int_interleaving_array)
properties = [];
N = length( int_interleaving_array );

n = length( int_pre_interleav );

if  N == n
    int_aft_interleav=zeros(n,1);
    for i = 1:N
        add = int_interleaving_array( i );
        int_aft_interleav( add ) = int_pre_interleav( i );
    end
else
    error('Receiving  wrong  array');
end