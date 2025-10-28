function [yri, ymi, y_fractional_outInt, ymi_HilbertInt, fractional_mult, fractional_sum, fractional_width_total_mult, ...
    fractional_width_total_sum, hilbert_mult, hilbert_sum, hilbert_width_total_mult, hilbert_width_total_sum] = ...
    fractional_delays(input_signal, hri_w, hh_m, coeff_frac_int, fractional_width, hilbert_coeff_int, enable_mask, fractional_mult_f, fractional_sum_f, ...
    hilbert_mult_f, hilbert_sum_f, sim_options)

        shift_frac_out = 13;
        shift_hilbert_out = 13;
        fractional_remainder = fractional_width - 1 - shift_frac_out;
        hilbert_remainder = 4;
        %%
        yri = filter(hri_w, 1, input_signal); % filter (стр.6 (15))
        ymi = filter(hh_m.', 1, yri);

        % y(n) = x(n)*(k1*2^N)+x(n-1)*(k2*2^N)+x(n-3)*(k3*2^N)
        % N = 18
        % y_fractional_outInt >> 13
        % N = 18 - 13 = 5;
        [y_fractional_outInt, fractional_mult, fractional_sum, fractional_width_total_mult, fractional_width_total_sum] = fir_filter(coeff_frac_int, input_signal, ...
            enable_mask, fractional_mult_f, fractional_sum_f, sim_options.width_fractional, sim_options); % (стр.6 (15)) 
        % y_fractional_outInt = bitshift(y_fractional_outInt,-shift_frac_out);
        % fractional_delay_filter_out_width = define_of_width_int(min(y_fractional_outInt)); % int18

        % if fractional_delay_filter_out_width > 18
        %     disp('Warning, data out fractional filter overflow!')
        %     disp([sim_options.SNR, sim_options.freq])
        % end

        snr_fractional_out_double = snr(yri, 1000000000);
        snr_fractional_out_int = snr(double(y_fractional_outInt)*2^-(18), 1000000000);

        if (snr_fractional_out_double - snr_fractional_out_int) > 0.1
            disp('SNR fractional out different!')
            disp([sim_options.SNR, sim_options.freq])
        end

        y_fractional_outInt_double = double(y_fractional_outInt)*2^-(18);
        relative_error_fractional = yri./y_fractional_outInt_double;


        figure(4);
        subplot(4,1,1)
        plot([yri(1:500), y_fractional_outInt_double(1:500)]);
        subplot(4,1,2);
        snr(yri, 1000000000);
        subplot(4,1,3);
        snr(y_fractional_outInt_double, 1000000000);


        subplot(4,1,4);
        plot(relative_error_fractional);
        title('Относительная ошибка выходного сигнала фильтра дробной задержки между double и integer')
        xlabel('Номер отсчета') 
        ylabel('Значение ошибки') 
        x4 = xline(37, '--', 'Переходной процесс фильтра')
        x4.LabelHorizontalAlignment = 'center'
        x4.LabelVerticalAlignment = 'middle';
        %% Hilbert
        %%
        % y(n) = x(n)*(k1*2^N)+x(n-1)*(k2*2^N)+x(n-3)*(k3*2^N)
        % N = 13
        [ymi_HilbertInt, hilbert_mult, hilbert_sum, hilbert_width_total_mult, hilbert_width_total_sum] = fir_filter(hilbert_coeff_int, y_fractional_outInt, ...
            enable_mask, hilbert_mult_f, hilbert_sum_f, sim_options.width_hilbert, sim_options); % (стр.6 (15)) );

        % ymi_HilbertInt = (bitshift(ymi_HilbertInt, -shift_hilbert_out));
        % hilbert_filter_out_width = define_of_width_int(min(ymi_HilbertInt)); % int18
        % % 
        % snr_hilbert_out_double = snr(ymi, 1000000000);
        % snr_hilbert_out_int = snr(double(ymi_HilbertInt)*2^-30, 1000000000);
        % 
        % 
        % if hilbert_filter_out_width > 18
        %     disp('Warning, data out Hilbert filter overflow!')
        %     disp([sim_options.SNR, sim_options.freq])
        % end
        % 
        % if (snr_hilbert_out_double - snr_hilbert_out_int) > 0.1
        %     disp('SNR hilbert out different!')
        %     disp([sim_options.SNR, sim_options.freq])
        % end
        % 
        figure(5);
        subplot(3,1,1)
        plot([ymi(1:500), double(ymi_HilbertInt(1:500))*2^-30]);
        subplot(3,1,2);
        snr(ymi, 1000000000);
        subplot(3,1,3);
        snr(double(ymi_HilbertInt)*2^-30, 1000000000);

end