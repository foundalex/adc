
function out = cordic_divide(numerator, denominator, iterations)


    scale = 2^11;

    xt = numerator;
    yt = denominator

    out = 0;

    for i = 0:iterations-1
        % против часовой стрелки
        if (xt <= 0)
            xn_t1 = xt + (fi(yt,1,70,55),-i,'int16'); % fi(1,70,55) + fi(1,70,55) = fi(1,71,55)
            out = out - bitshift(scale,-i,'int16');
        else
            % по часовой стрелке
            xn_t1 = xt - bitshift(int16(yt),-i,'int16'); % fi(1,13,10) + fi(1,12,10) = fi(1,14,10)
            out = out + bitshift(scale,-i,'int16');
        end

        % xn_t1 = fi(bitshift(int16(xn_t1),-2,'int16'),1,12,0); % fi(1,12,8) 
        xt = xn_t1;

    end

end