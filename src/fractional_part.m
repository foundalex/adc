function [yri, yri_round, ymi, ymi_round, ...
    y_fractional_outInt_div, filter_max_width_out_fractional, ...
    ymi_HilbertInt_div, filter_max_width_out_hilbert] = fractional_part ...
    (   bandpass_fractional, bandpass_hilbert, ...
        coeff_frac_int, hilbert_coeff_int, ...
        y_golden_outInt_div, ...
        yr_double, del_proc, fractional_max_width, hilbert_max_width, sim_options)


        filter_max_width_out_fractional = struct;
        filter_max_width_out_hilbert = struct;


        %% Фильтр дробной задержки (double)
        yri = filter(bandpass_fractional, 1, yr_double); % filter (стр.6 (15))
        % убираем переходной процесс
        yri = [yri(del_proc+1:end); zeros(del_proc,1)]; % убираем переходной процесс
        % округляем значения
        yri_round = round(yri);

        %% Фильтр Гилберта (Double)
        ymi = filter(bandpass_hilbert.', 1, yri);
        % убираем переходной процесс
        ymi = [ymi(del_proc+1:end); zeros(del_proc,1)]; 
        % округляем значения
        ymi_round = round(ymi);
        
        %% Фильтр дробной задержки (int)
        % [y_fractional_outInt, filter_max_width_out_fractional] = ...
        %     fir_filter(coeff_frac_int, y_golden_outInt_div, sim_options.N, fractional_max_width, sim_options.width_fractional, sim_options); % (стр.6 (15)) 
     
        y_fractional_outInt = int64(filter(coeff_frac_int, int64(1), y_golden_outInt_div)); % filter (стр.6 (15))


        % убираем переходной процесс
        y_fractional_outInt = [y_fractional_outInt(del_proc+1:end); zeros(del_proc,1)]; 
        % округляем значения
        y_fractional_outInt_div = round_int(y_fractional_outInt, sim_options.fractional_coeff_width-1, sim_options.type_fir_out);


        %% Фильтр Гилберта (int)
        % [ymi_HilbertInt, filter_max_width_out_hilbert] = ...
        %     fir_filter(hilbert_coeff_int, y_fractional_outInt_div, sim_options.N, hilbert_max_width, sim_options.width_hilbert, sim_options); % (стр.6 (15)) );
        % ymi_HilbertInt = [ymi_HilbertInt(del_proc+1:end); zeros(del_proc,1)]; % убираем переходной процесс

        ymi_HilbertInt = int64(filter(hilbert_coeff_int, int64(1), y_fractional_outInt_div)); % filter (стр.6 (15))

        % убираем переходной процесс
        ymi_HilbertInt = [ymi_HilbertInt(del_proc+1:end); zeros(del_proc,1)]; 

        % округляем значения
        ymi_HilbertInt_div = round_int(ymi_HilbertInt, sim_options.hilbert_coeff_width-1, sim_options.type_fir_out);

end