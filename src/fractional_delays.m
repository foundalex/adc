function [yri, ymi, y_fractional_outInt, ymi_HilbertInt, hilbert_out_width] = fractional_delays(input_signal, Z, hri_w, hh_m, fractional_width, hilbert_width)

        frac_out_width = 32;
        hilbert_out_width = 63;
        %%
        yri = filter(hri_w, 1, input_signal(:,1)); % filter (стр.6 (15))
        ymi = filter(hh_m.', 1, yri);

        %% integer
        coeff_frac_int = fi(hri_w*2^(fractional_width-1),1,fractional_width,0);

       
        y_fractional_outInt = fir_filter(coeff_frac_int, input_signal(:,1), 12, frac_out_width); % (стр.6 (15)) 
        if (double(y_fractional_outInt) == 2^32)
            exit;
        end
        % fi(1,19,18) * fi(1,12,0) = fi(1,31,18) + 1 = fi(1,32,18)

        % y_fractional_outInt = bitshift(y_fractional_outInt,-15);
        % y_fractional_outInt = fi(y_fractional_outInt,1,17,0);
        % fi(1,32,18) - 15 = 1,17,3

        % figure(4);
        % subplot(3,1,1)
        % plot([yri(1:250), double(y_fractional_outInt(1:250))*2^-(3)]);
        % subplot(3,1,2);
        % snr(yri, 8000000000);
        % subplot(3,1,3);
        % snr(double(y_fractional_outInt)*2^-(3), 8000000000);

        %% Hilbert
        % Negative Symmetric coefficients
        hh_m_int = fi(hh_m*2^(hilbert_width-1),1,hilbert_width,0); 

        for i = 1:73
            for j = 1:64
                if (abs(hh_m_int(i)) < 2^j)
                    if (hh_m_int(i)) < 0
                        width(i) = j+1;
                        break;
                    % elseif (hh_m_int(i) == 0)
                    %     width(i) = 0;
                    %     break;
                    else
                        width(i) = j;
                        break;
                    end
                elseif (abs(hh_m_int(i)) == 2^j)
                    width(i) = j+1;
                end
            end
        end

        figure(3)
        subplot(2,1,1)
        plot(double(hh_m_int))
        title('Импульсная характеристика фильтра Гилберта')
        xlabel('Номер отсчета') 
        ylabel('Амплитуда') 
        subplot(2,1,2)
        plot(width)
        title('Разрядность коэффициентов')
        xlabel('Номер коэффициента') 
        ylabel('Необходимое количество бит') 

        ymi_HilbertInt = fir_filter(hh_m_int, y_fractional_outInt, frac_out_width, hilbert_out_width); % (стр.6 (15)) 
        % fi(1,13,12) * fi(1,32,18) = fi(1,45,30) + 18 = fi(1,63,30)
        % ymi_HilbertInt = (bitshift(ymi_HilbertInt, -12));
        % ymi_HilbertInt = fi(ymi_HilbertInt,1,36,0);
        % % fi(1,13,12) * fi(1,17,3) = fi(1,30,15) + 18 = fi(1,48,15)
        % % fi(1,63,30) - 12 = fi(1,51,18)

        % figure(5);
        % subplot(3,1,1)
        % plot([ymi(1:250), double(ymi_HilbertInt(1:250))*2^-30]);
        % subplot(3,1,2);
        % snr(ymi, 8000000000);
        % subplot(3,1,3);
        % snr(double(ymi_HilbertInt)*2^-30, 8000000000);

end