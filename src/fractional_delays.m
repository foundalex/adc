function [yri, ymi, y_fractional_outInt, ymi_HilbertInt] = fractional_delays(input_signal, Z, hri_w, hh_m, fractional_width, hilbert_width)

       
        yri = filter(hri_w, 1, input_signal(:,1)); % filter (стр.6 (15))

        %% integer
        hrim_int = fi(hri_w*2^(fractional_width-1),1,fractional_width,0);

        y_fractional_outInt = fir_filter(hrim_int, input_signal(:,1), 12, 33); % (стр.6 (15)) 
        if (double(y_fractional_outInt) == 2^32)
            exit;
        end
        % fi(1,19,18) * fi(1,12,0) = fi(1,31,18) + 1 = fi(1,33,18)


        % figure(4);
        % subplot(3,1,1)
        % plot([yri(1:250), double(y_fractional_outInt(1:250))*2^-(fractional_width-1)]);
        % subplot(3,1,2);
        % snr(yri, 8000000000);
        % subplot(3,1,3);
        % snr(double(y_fractional_outInt)*2^-(fractional_width-1), 8000000000);


        ymi = filter(hh_m.', 1, yri);

        hh_m_int = fi(hh_m*2^(hilbert_width-1),1,hilbert_width,0);

        % hh_m_int_mix = zeros(73,1);
        % hh_m_int_mix(1:2:73) = hh_m_int(1:37);
        % hh_m_int_mix(2:2:72) = hh_m_int(38:73);

        ymi_HilbertInt = fir_filter(hh_m_int, y_fractional_outInt, 33, 62); % (стр.6 (15)) 

        maxHilbert = max(ymi_HilbertInt);
        % fi(1,13,12) * fi(1,33,18) = fi(1,46,30) + 16 = fi(1,62,30)

        % ymi_HilbertInt = (bitshift(ymi_HilbertInt, -12)); % fi(1,19,12)
        % 
        % figure(5);
        % subplot(3,1,1)
        % plot([ymi(1:250), double(ymi_HilbertInt(1:250))*2^-30]);
        % subplot(3,1,2);
        % snr(ymi, 8000000000);
        % subplot(3,1,3);
        % snr(double(ymi_HilbertInt)*2^-30, 8000000000);

end