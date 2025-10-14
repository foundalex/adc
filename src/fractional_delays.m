function [yri, ymi, y_fractional_outInt, ymi_HilbertInt] = fractional_delays(input_signal, Z, hri_w, hh_m, fractional_width, hilbert_width)

        frac_out_width = 31;
        hilbert_out_width = 30;
        %%
        yri = filter(hri_w, 1, input_signal(:,1)); % filter (стр.6 (15))
        ymi = filter(hh_m.', 1, yri);

        %% integer
        coeff_frac_int = fi(hri_w*2^(fractional_width-1),1,fractional_width,0);


        % coeff_frac_int(2) = fi(0,1,19,0);
        % coeff_frac_int(3) = fi(0,1,19,0);
        % coeff_frac_int(4) = fi(1,1,19,0);
        % 
        % a = input_signal(:,1)* fi(2,1,19,0);
        % a1 = bitshift(a,-4);
        % a1 = fi(a1,1,12,0);
        % 
        % b = a1* fi(-6,1,19,0);
        % a2 = bitshift(b,-4);
        % a2 = fi(a1,1,12,0);
        % 
        % c = a2* fi(13,1,19,0);
        % a3 = bitshift(c,-4);
        % a3 = fi(a1,1,12,0);
        % 
        % e = a1 + a2 + a3;


        % for i = 1:73
        %     for j = 1:64
        %         if (abs(coeff_frac_int(i)) < 2^j)
        %             if (coeff_frac_int(i)) < 0
        %                 width(i) = j+1;
        %                 break;
        %             % elseif (hh_m_int(i) == 0)
        %             %     width(i) = 0;
        %             %     break;
        %             else
        %                 width(i) = j;
        %                 break;
        %             end
        %         elseif (abs(coeff_frac_int(i)) == 2^j)
        %             width(i) = j+1;
        %         end
        %     end
        % end
        % 
        % figure(3)
        % subplot(2,1,1)
        % plot(double(coeff_frac_int))
        % title('Импульсная характеристика фильтра Гилберта')
        % xlabel('Номер отсчета') 
        % ylabel('Амплитуда') 
        % subplot(2,1,2)
        % plot(width)
        % title('Разрядность коэффициентов')
        % xlabel('Номер коэффициента') 
        % ylabel('Необходимое количество бит') 


        [y_fractional_outInt, y_fractional_outInt11] = fir_filter(coeff_frac_int, input_signal(:,1), 12, frac_out_width); % (стр.6 (15)) 
        % fi(1,19,18) * fi(1,12,0) = fi(1,31,18)

        % y_fractional_outInt = bitshift(y_fractional_outInt,-15);
        % y_fractional_outInt = fi(y_fractional_outInt,1,16,0);
        % fi(1,31,18) - 15 = 1,16,3

        figure(4);
        subplot(3,1,1)
        plot([yri(1:250), double(y_fractional_outInt11(1:250))']);
        subplot(3,1,2);
        snr(yri, 8000000000);
        subplot(3,1,3);
        snr(double(y_fractional_outInt11), 8000000000);

        %% Hilbert
        % Negative Symmetric coefficients
        hh_m_int = fi(hh_m*2^(hilbert_width-1),1,hilbert_width,0); 

        % for i = 1:73
        %     for j = 1:64
        %         if (abs(hh_m_int(i)) < 2^j)
        %             if (hh_m_int(i)) < 0
        %                 width(i) = j+1;
        %                 break;
        %             % elseif (hh_m_int(i) == 0)
        %             %     width(i) = 0;
        %             %     break;
        %             else
        %                 width(i) = j;
        %                 break;
        %             end
        %         elseif (abs(hh_m_int(i)) == 2^j)
        %             width(i) = j+1;
        %         end
        %     end
        % end
        % 
        % figure(3)
        % subplot(2,1,1)
        % plot(double(hh_m_int))
        % title('Импульсная характеристика фильтра Гилберта')
        % xlabel('Номер отсчета') 
        % ylabel('Амплитуда') 
        % subplot(2,1,2)
        % plot(width)
        % title('Разрядность коэффициентов')
        % xlabel('Номер коэффициента') 
        % ylabel('Необходимое количество бит') 

        ymi_HilbertInt = fir_filter(hh_m_int, y_fractional_outInt, frac_out_width, hilbert_out_width); % (стр.6 (15)) 
        % fi(1,13,12) * fi(1,16,3) = fi(1,29,15)
        ymi_HilbertInt = (bitshift(ymi_HilbertInt, -13));
        ymi_HilbertInt = fi(ymi_HilbertInt,1,16,0);
        % % fi(1,29,15) - 13 = fi(1,16,2)

        figure(5);
        subplot(3,1,1)
        plot([ymi(1:250), double(ymi_HilbertInt(1:250))*2^-2]);
        subplot(3,1,2);
        snr(ymi, 8000000000);
        subplot(3,1,3);
        snr(double(ymi_HilbertInt)*2^-2, 8000000000);

end