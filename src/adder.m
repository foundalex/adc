function [y, overflow, width_total] = adder(a,b,int_size,N)

    %% проверка выходной разрядности
    if (a < 0)
        a_abs = a * cast(-1, int_size);
    else
        a_abs = a;
    end

    if (b < 0)
        b_abs = b * cast(-1, int_size);
    else
        b_abs = b;
    end

    width_a = define_of_width_int(a_abs,int_size);
    width_b = define_of_width_int(b_abs,int_size);

    width_total = 0; %width_a + width_b - 1;

    %% сумматор

    y = a + b;

    if (y > (2^(N-1)-1) | y < -2^(N-1))
        overflow = int8(1);
    else
        overflow = int8(0);
    end


end