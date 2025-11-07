function width = define_of_width_int(input, int_size, num)
    % находим модуль числа
    if (input < 0)
        input_abs = input * cast(-1, int_size);
    else
        input_abs = input;
    end

    if int_size == "int64"
        for i = 0:num-1
            table_list(i+1) = cast(2^i,int_size);
        end
        %%
        if (input_abs == 0)
            width = cast(0, int_size);
        else
            for i = cast(1:num,int_size)
                mask = table_list(num-i);
                out1 = bitand(mask,input_abs);
                if (out1 > 0)
                    break;
                end
            end
            width = num-i+cast(1,int_size); % добавляем 1 разряд для знака
        end
    elseif int_size == "double" || int_size == "single" 
        for i = 1:90
            if ((2^i)-1 >= input_abs)
                width = i+1; % 1 бит для знака
                break;
            end
        end
    end
end