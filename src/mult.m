function [y, overflow, width_total] = mult(a, b, int_size, width_fractonal)
    % проверка выходной разрядности

    width_a = define_of_width_int(a, int_size, width_fractonal);
    width_b = define_of_width_int(b, int_size, width_fractonal);

    width_total = width_a + width_b - 1;

    %% умножитель
    y = a * b;

    if (y >= (2^(width_fractonal-1))-1 | y < -2^(width_fractonal-1))
        overflow = cast(1,int_size);
    else
        overflow = cast(0,int_size);
    end

end