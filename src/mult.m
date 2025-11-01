function [y, overflow, y_abs, width_total] = mult(a, b, int_size, width)
    % проверка выходной разрядности

    % if int_size == "int64"
        width_a = define_of_width_int(a, int_size, width);
        width_b = define_of_width_int(b, int_size, width);

        width_total = width_a + width_b - 1;
    % else
    %     width_total = 80;
    % end

    %% умножитель
    mmax = cast((2^(width-1))-1, int_size);
    mmin = cast(-2^(width-1), int_size);

    y = a * b;

    if (y >= mmax || y <= mmin)
        overflow = cast(1,int_size);
    else
        overflow = cast(0,int_size);
    end

    %% число по модулю
    if y < 0
        y_abs = y * cast(-1, int_size); % находим число по модулю
    else
	    y_abs = y;
    end
end