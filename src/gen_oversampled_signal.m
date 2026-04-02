function [s_to_subadc, adc_input_int, s_after_subadc] = gen_oversampled_signal(sim_options)

    dt = 1/sim_options.Fs;                                                              % seconds per sample
    t = 0:dt:sim_options.StopTime;                                                      % seconds

    % Create main signal with noise in double
    s1 = sin(2*pi*sim_options.freq*t);

    s = awgn(s1,sim_options.SNR(1), 'measured');

    % шаг выборки из передискретизированного сигнала
    step = sim_options.M*sim_options.Inter;
    % начало выборки
    begin = (1:sim_options.M).*sim_options.Inter-(sim_options.Inter-1);

    if sim_options.MODEL_ERROR == true
        dlin = floor((length(s)-(sim_options.time_skew_array(1:sim_options.M-1).*sim_options.Inter))./(sim_options.M*sim_options.Inter));
        dlin = dlin - 1;
    else
        dlin = (floor((length(s)-begin)./(sim_options.M*sim_options.Inter)));
    end

    ended = max(dlin)*sim_options.Inter*sim_options.M;
    
    %% разбиваем входной сигнал на сигналы для суб-АЦП
    
    for i = 1:sim_options.M  
        adc_input_good(:,i) = s(begin(i):step:ended);
        adc_input(:,i) = s(begin(i):step:ended);

        % добавление ошибок к каналам АЦП
        if sim_options.MODEL_ERROR == true
            % time skew
            adc_input(:,i) = s(begin(i) + sim_options.time_skew_array(i)*sim_options.Inter:step:ended);
            % offset
            adc_input(:,i) = adc_input(:,i) + sim_options.offset_error_array(i);
            % gain
            adc_input(:,i) = adc_input(:,i) * sim_options.gain_error_array(i);
        end
    end

    % перевод в инты
    for i = 1:sim_options.M
        adc_input_int(:,i) = int16(round(fi(adc_input(:,i),1,12,11)*sim_options.Bit));
        adc_input_int_good(:,i) = int16(round(fi(adc_input_good(:,i),1,12,11)*sim_options.Bit));
    end

	% итоговый выходной сигнал до искажений
	s_to_subadc = zeros(sim_options.M*length(adc_input_int_good(:,1)),1);
    % итоговый выходной сигнал после искажений
    s_after_subadc = zeros(sim_options.M*length(adc_input(:,1)),1);

	for i = 1:sim_options.M
		s_to_subadc(i:sim_options.M:end) = adc_input_int_good(:,i); 
        s_after_subadc(i:sim_options.M:end) = adc_input_int(:,i); 
    end

    % ЧХ выходных сигналов
    figure(3);
    subplot(4,1,1)
    plot(s_to_subadc(1:100));
    subplot(4,1,2)
    plot(s_after_subadc(1:100));
    subplot(4,1,3)
    snr(s_to_subadc, sim_options.Fs/sim_options.Inter);
    subplot(4,1,4)
    snr(s_after_subadc, sim_options.Fs/sim_options.Inter);

    % Выходные отсчеты каналов АЦП
    figure(54);
    subplot(sim_options.M+1,1,1)
    plot(s_after_subadc(1:100));
    title('Выходной сигнал всего АЦП');
    xlabel('Амплитуда');
    ylabel('Номер отсчета'); 
    for i = 1:sim_options.M
        subplot(sim_options.M+1,1,i+1)
        plot(adc_input_int(1:100,i));
        title(['Выход АЦП' num2str(i)])
        xlabel('Амплитуда') 
        ylabel('Номер отсчета') 
    end

end