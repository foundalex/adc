function out = bitmask(input, width, int_size)

    if input < 0
        m = cast(-1, int_size);
    else
        m = cast(1, int_size);
    end

    mask = cast((2^width)-1,int_size);
    abs_input = input * m;

    out = bitand(mask,abs_input);

    out = out * m;

    % ww = int_size - width;
    % c = bitshift(input, ww);
    % 
    % out = c / (2^ww);


end