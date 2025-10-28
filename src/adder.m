function [y, overflow, y_abs, width_total] = adder(a,b,int_size,width)

    % проверка выходной разрядности

    width_a = define_of_width_int(a,int_size,width);
    width_b = define_of_width_int(b,int_size,width);

    if (a <= 0 && b <= 0) || (a > 0 && b > 0) % если операнды с одинаковыми знаками, то прибавляем к самому большому 1
        if width_a >= width_b
            width_total = width_a + cast(1,int_size); 
        else
            width_total = width_b + cast(1,int_size);
        end
    else
        if width_a >= width_b % если с разными  
            width_total = width_a; 
        else
            width_total = width_b;
        end
    end

    %% сумматор
    mmax = cast((2^(width-1))-1, int_size);
    mmin = cast(-2^(width-1), int_size);

    y = a + b;

    if (y > mmax | y < mmin)
        overflow = int8(1);
    else
        overflow = int8(0);
    end

    %% число по модулю
    if y < 0
        y_abs = y * cast(-1, int_size); % находим число по модулю
    else
        y_abs = y;
    end

end