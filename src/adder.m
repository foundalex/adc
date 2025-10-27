function [y, overflow, width_total] = adder(a,b,int_size,width_fractonal)

    % проверка выходной разрядности

    width_a = define_of_width_int(a,int_size,width_fractonal);
    width_b = define_of_width_int(b,int_size,width_fractonal);

    if (a <= 0 && b <= 0) || (a > 0 && b > 0) % если операнды с одинаковыми знаками, то прибавляем к самому большому 1
        if width_a >= width_b
            width_total = width_a + 1; 
        else
            width_total = width_b + 1;
        end
    else
        if width_a >= width_b % если с разными  
            width_total = width_a; 
        else
            width_total = width_b;
        end
    end

    %% сумматор
    y = a + b;

    if (y > cast((2^(width_fractonal-1)-1),int_size) | y < cast(-2^(width_fractonal-1),int_size))
        overflow = int8(1);
    else
        overflow = int8(0);
    end

end