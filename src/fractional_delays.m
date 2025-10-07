function [yri, ymi, y_fractional_outInt, ymi_HilbertInt] = fractional_delays(input_signal, input_signal_int, Z, hri_w, hh_m, fractional_width, hilbert_width)

       
        yri = filter(hri_w, 1, input_signal(:,1)); % filter (стр.6 (15))

        %% integer
        hrim_fi = fi(hri_w, 1,fractional_width,fractional_width-1);
        hrim_int = int32(hrim_fi * 2^(fractional_width-1));

        y_fractional_outInt = fir_filter(hrim_int, input_signal_int(:,1)); % (стр.6 (15)) 
        % y_fractional_outInt = int32(bitshift(y_fractional_outInt, 0)); % fi(1,13,11)

        % figure(4);
        % plot([yri(1:250), double(y_fractional_outInt16(1:250))*2^0]);
        % % 
        % figure(10);
        % subplot(2,1,1);
        % snr(yri, 8000000000);
        % subplot(2,1,2);
        % snr(double(y_fractional_outInt16)*2^0, 8000000000);


        ymi = filter(hh_m.', 1, yri);

        hh_m_fi = fi(hh_m,1,hilbert_width,hilbert_width-1);
        hh_m_int = int16(round(hh_m_fi * 2^(hilbert_width-1))); % int coeff fi(1,4,3)

        ymi_HilbertInt = fir_filter(hh_m_int, y_fractional_outInt); % (стр.6 (15)) fi(1,16,15) * fi(1,14,11) = fi(1,30,27)
        % ymi_int16 = int16(round(ymi_HilbertInt/32768)); % fi(1,30,27) - 16 bit = fi(1,14,11)
        ymi_HilbertInt = ymi_HilbertInt; % (bitshift(ymi_HilbertInt, -1)); % fi(1,13,11)


        % figure(4);
        % plot([ymi(1:250), double(ymi_HilbertInt(1:250))*2^-12]);
        % % 
        % figure(10);
        % subplot(2,1,1);
        % snr(ymi, 8000000000);
        % subplot(2,1,2);
        % snr(double(ymi_HilbertInt)*2^-12, 8000000000);

        % y_out_Hilbert_max(:,i) = max(ymi_HilbertInt);
        % y_out_Hilbert_min(:,i) = min(ymi_HilbertInt);

        % max_fractional = max(y_fractional_outInt_max);
        % min_fractional = min(y_fractional_outInt_min);

end