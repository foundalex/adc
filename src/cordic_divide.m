
function out = cordic_divide(numerator, denominator, iterations)

    scale = 2^11;

    xt = numerator;
    yt = denominator;

    out = 0;

    for i = 0:iterations-1
        % против часовой стрелки
        if (xt <= 0)
            xn_t1 = int16(xt) + bitshift(yt,-i,'int16');
            out = out - bitshift(scale,-i,'int16');
        else
            % по часовой стрелке
            xn_t1 = int16(xt) - bitshift(yt,-i,'int16');
            out = out + bitshift(scale,-i,'int16');
        end
        xt = xn_t1;
    end

end