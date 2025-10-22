function [c, a] = bitmask(input, width, int_size)

    if input < 0
        m = cast(-1, int_size);
    else
        m = cast(1, int_size);
    end

    a = input * m;
    c = bitand(cast(2^width-1,int_size),a);


    % a = bitget(int32(bitmask_i), width);
    % if a == 1
    %     c = (int32(bitmask_i)) - int32(2^width);
    % else
    %     c = input;
    % end


end