function [int_output] = de_interleaving( int_channeloutput_interleav ,int_interleaving_array )

N = length( int_interleaving_array );

n = length( int_channeloutput_interleav );

if  N == n
	int_output=zeros(n,1);   
	for i = 1:N
	    add = int_interleaving_array( i );
	    int_output( i ) = int_channeloutput_interleav( add );
    end    
else 
   error('receiving wrong array');
end