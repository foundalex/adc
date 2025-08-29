
function out = cordic_divide(numerator, denominator, iterations)


    scale = 2^11;

    xt = fi(numerator,1,13,0); % fi(1,13,10)
    yt = fi(denominator,1,12,0); % fi(1,12,10)

    out = 0;

    for i = 0:iterations-1
        % против часовой стрелки
        if (xt <= 0)
            xn_t1 = xt + fi(bitshift(int16(yt),-i,'int16'),1,12,0); % fi(1,13,10) + fi(1,12,10) = fi(1,14,10)
            out = out - bitshift(scale,-i,'int16');
        else
            % по часовой стрелке
            xn_t1 = xt - fi(bitshift(int16(yt),-i,'int16'),1,12,0); % fi(1,13,10) + fi(1,12,10) = fi(1,14,10)
            out = out + bitshift(scale,-i,'int16');
        end

        % xn_t1 = fi(bitshift(int16(xn_t1),-2,'int16'),1,12,0); % fi(1,12,8) 
        xt = xn_t1;

    end

end