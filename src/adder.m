function [y, overflow] = adder(a,b,N)


    y = a + b;

    if (y > (2^N)-1)
        overflow = 1;
    else
        overflow = 0;
    end


end