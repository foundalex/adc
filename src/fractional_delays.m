function [yri, ymi, y_fractional_outInt, ymi_HilbertInt, fractional_mult, fractional_sum, hilbert_mult, hilbert_sum] = fractional_delays(input_signal, hri_w, hh_m, coeff_frac_int, ...
    fractional_width, hilbert_coeff_int, enable_mask, fractional_mult_f, fractional_sum_f, sim_options)

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
        [y_fractional_outInt, fractional_mult, fractional_sum] = fir_filter(coeff_frac_int, input_signal, 'int32', int32(32), ...
            enable_mask, fractional_mult_f, fractional_sum_f, sim_options); % (стр.6 (15)) 
        % y_fractional_outInt = bitshift(y_fractional_outInt,-shift_frac_out);
        fractional_delay_filter_out_width = define_of_width_int(min(y_fractional_outInt)); % int18

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

        figure(4);
        subplot(3,1,1)
        plot([yri(1:250), double(y_fractional_outInt(1:250))*2^-(18)]);
        subplot(3,1,2);
        snr(yri, 1000000000);
        subplot(3,1,3);
        snr(double(y_fractional_outInt)*2^-(18), 1000000000);

        %% Hilbert
        %%
        ymi_HilbertInt = 0;
        hilbert_mult = 0;
        hilbert_sum = 0;
        % y(n) = x(n)*(k1*2^N)+x(n-1)*(k2*2^N)+x(n-3)*(k3*2^N)
        % N = 13
        % y_fractional_outInt = 1,18,5 * 1,13,12 = 1,31,17;
        % y_fractional_outInt >> 13
        % N = 1,31,17 - 13 = 1,18,4;
        % [ymi_HilbertInt, hilbert_mult, hilbert_sum] = fir_filter(int64(hilbert_coeff_int), y_fractional_outInt, 'int64', 64, ...
        %     0, 'Width multiplier Hilbert filter.txt', 'Width adder Hilbert filter.txt'); % (стр.6 (15)) ); % (стр.6 (15)) 
        % % ymi_HilbertInt = (bitshift(ymi_HilbertInt, -shift_hilbert_out));
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
        % figure(5);
        % subplot(3,1,1)
        % plot([ymi(1:250), double(ymi_HilbertInt(1:250))*2^-30]);
        % subplot(3,1,2);
        % snr(ymi, 1000000000);
        % subplot(3,1,3);
        % snr(double(ymi_HilbertInt)*2^-30, 1000000000);

end