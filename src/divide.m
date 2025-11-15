function [efi_out, overflow, c_abs, width_total] = divide(a, b, int_size, width, factor)

    % проверка выходной разрядности
    width_a = define_of_width_int(a, int_size, width);
    width_b = define_of_width_int(b, int_size, width);

    if (width_a > width_b)
        width_total = width_a;
    else
        width_total = width_b;
    end

    %% делитель
    mmax = cast((2^(width-1))-1, int_size);
    mmin = cast(-2^(width-1), int_size);

    a1 = double(a);
    b1 = double(b);
    c = a1 / b1;
  
    % efi_out = double(c*2^width);
    efi_out = c;

    if (c >= mmax || c <= mmin)
        overflow = cast(1,int_size);
    else
        overflow = cast(0,int_size);
    end
    %% число по модулю
    if c < 0
        c_abs = c * cast(-1, int_size); % находим число по модулю
    else
	    c_abs = c;
    end
end
