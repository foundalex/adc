function out = bitmask(input, int_size, width)
    % 

    table_list = cast(zeros(width+1,1),int_size);

    for i = 1:width+1
        table_list(i) = cast((2^(i+1))-1,int_size);
    end

    %%
    if input < 0
        m = cast(-1, int_size);
    else
        m = cast(1, int_size);
    end
    %%
    if width == 0
        out = 0;
    else
        mask = table_list(width-2);
        abs_input = input * m;
        % 
        out1 = bitand(mask,abs_input);
    
        if (out1 == 0)
            out = input;
        else
            out = out1 * m;
        end
    end

end

 