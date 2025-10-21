function [y, overflow] = mult(a,b,N)

    y = a * b;

    if (y >= (2^(N-1))-1 | y < -2^(N-1))
        overflow = 1;
        disp('overflow mult detected');
    else
        overflow = 0;
    end

end