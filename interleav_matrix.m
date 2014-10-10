%% generate a pseudo-random interleaving matrix

function [int_interleaving_array] = interleav_matrix( int_interleaving_buffer )

N = length(int_interleaving_buffer);
int_interleaving_array = [];

while ( length( int_interleaving_array ) ~= N )
    
    int_interleaving_array = [];
    
    int_number_array = 1:N;

    int_refuse_array = []; 
    
    while  length( int_number_array ) > 0
        
        n = length(int_number_array);
        
        subscript = floor( n * rand + 1 );

        int_rand_num = int_number_array( subscript );

        int_number_array( subscript ) = int_number_array( n );

        if  n == 1
            int_number_array = [];
        else
            int_number_array = int_number_array( 1:n-1 );
        end
        
        s = floor( (N/2)^0.5 ) - 1;
        flag = 1;

        len = length(int_interleaving_array);
        if  len < s
            s = len;
        end

        if  n <= 2 * s
            s = floor( n/2 );
        end

        for i = 0:s-1
            if abs( int_interleaving_array ( len - i ) - int_rand_num ) <= s
                flag = 0;
                break;
            end
        end
        
        if  flag == 0
            int_refuse_array = int_rand_num;
        else
            int_interleaving_array = [int_interleaving_array,int_rand_num];
        end
        
        if  length( int_refuse_array ) > 0
            int_number_array = [int_number_array,int_refuse_array];
            int_refuse_array = [];
        end
    end
end
            
        