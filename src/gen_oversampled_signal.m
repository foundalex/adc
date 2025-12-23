function [s_to_subadc, adc_input_double, adc_input_int, s_after_subadc, Z] = gen_oversampled_signal(sim_options)

    dt = 1/sim_options.Fs;                                                              % seconds per sample
    t = 0:dt:sim_options.StopTime;                                                      % seconds

    % Create main signal with noise in double
    
    s1 = (cos(2*pi*sim_options.freq*t));
    % s2 = cos(2*pi*(sim_options.freq+10000000)*t+pi/4);
    % s3 = (cos(2*pi*(sim_options.freq+20000000)*t));
    % s4 = (cos(2*pi*(sim_options.freq+30000000)*t));
    % s5 = (cos(2*pi*(sim_options.freq+40000000)*t));

    % figure(3);
    % plot([s1(1:100)', s2(1:100)']);

    s = s1; %  + s2; % + s3 + s4 + s5;
    % s = awgn(s, 60);
    % s = s + noise;
    % [pxx,f] = periodogram(s); 

    step = sim_options.M*sim_options.Inter;
    begin = (1:sim_options.M).*sim_options.Inter+1;

    if sim_options.MODEL_ERROR == true
        dlin = floor((length(s)-(sim_options.time_skew_array(1:sim_options.M-1).*sim_options.Inter))./(sim_options.M*sim_options.Inter));
        dlin = dlin - 1;
    else
        dlin = floor((length(s)-begin)./(sim_options.M*sim_options.Inter));
    end

    adc_input = zeros(dlin(sim_options.M-1),sim_options.M);
    adc_input_int = int16(zeros(dlin(sim_options.M-1),sim_options.M));

    ended = dlin(sim_options.M-1)*sim_options.Inter*sim_options.M;
    
    %% разбиваем входной сигнал на сигналы для суб-АЦП

    % АЦП0 - эталон
    adc_input(:,1) = s(begin(1):step:ended);
    adc_input_good(:,1) = s(begin(1):step:ended);
    
    for i = 2:sim_options.M  
        % если ошибки суб-АЦП включены, то добавляем time skew
        if i == sim_options.M
            adc_input_good(:,i) = s(begin(i):step:ended+sim_options.Inter);
            if sim_options.MODEL_ERROR == true
                adc_input(:,i) = s(begin(i) + sim_options.time_skew_array(i-1)*sim_options.Inter:step:ended+10);
            else
                adc_input(:,i) = s(begin(i):step:end-1);
            end
        else
            adc_input_good(:,i) = s(begin(i):step:ended);
            if sim_options.MODEL_ERROR == true
                adc_input(:,i) = s(begin(i) + sim_options.time_skew_array(i-1)*sim_options.Inter:step:ended+10);
            else
                adc_input(:,i) = s(begin(i):step:ended);
            end
        end
    end

    %% добавляем к каждому суб-АЦП шум
    adc_input_double = adc_input;
    for i = 1:sim_options.M 
        adc_input(:,i) = awgn(adc_input(:,i), sim_options.SNR(i) , "measured");
        adc_input_double(:,i) = awgn(adc_input_double(:,i), sim_options.SNR(i) , "measured");
        adc_input_good(:,i) = awgn(adc_input_good(:,i), sim_options.SNR(i) , "measured");
    end

    for i = 1:sim_options.M 
        s_fi(:,i) = fi(adc_input(:,i),1,12,11);
        s_int(:,i) = int16(round(s_fi(:,i)*2^11));
    end

    %% если ошибки суб-АЦП включены, то добавляем gain error
    adc_input_int(:,1) = s_int(:,1);
    % adc_input_double(:,1) = adc_input(:,1);
    if sim_options.MODEL_ERROR == true
        for i = 1:sim_options.M-1 
            adc_input_double(:,i+1) = (adc_input(:,i+1) * sim_options.gain_error_array(i));
            adc_input_int(:,i+1) = int16(fi(double(s_int(:,i+1)) * sim_options.gain_error_array(i),1,12,0));
        end
    end

    % % model offset error
    % for i = 1:sim_options.M-1
    %     adc_input(:,i+1) = adc_input(:,i+1) + offset_error_array(i);
    % end

    % перевод в инты
    for i = 1:sim_options.M
        adc_input_int(:,i) = int16(round(fi(adc_input(:,i),1,12,11)*sim_options.Bit));
        adc_input_int_good(:,i) = int16(round(fi(adc_input_good(:,i),1,12,11)*sim_options.Bit));
        % adc_input_int_good(:,i) = bitshift(adc_input_int_good(:,i), -2);
    end


    % spectrumScope = spectrumAnalyzer(SampleRate=sim_options.Fs, ...            
    %         AveragingMethod='exponential',ForgettingFactor=0.99, ...
    %         YLimits=[-30 10],ShowLegend=true);
    % 
    % spectrumScope([double(adc_input_int(:,1))/2048]);
    % % spectrumScope([s']);


	% исходный сигнал до искажений
	s_to_subadc = zeros(sim_options.M*length(adc_input_int_good(:,1)),1);
	for i = 1:sim_options.M
		s_to_subadc(i:sim_options.M:end) = adc_input_int_good(:,i); 
	end

    if sim_options.MODEL_ERROR == false
        adc_input_int = adc_input_int_good;
    end

    % сигнал после искажений
    s_after_subadc = zeros(sim_options.M*length(adc_input_int(:,1)),1);
	for i = 1:sim_options.M
        s_after_subadc(i:sim_options.M:end) = adc_input_int(:,i); 
    end

    Z = ceil(sim_options.freq/(sim_options.Fs/sim_options.Inter/2/sim_options.M));      % Nyquist zone

    % save (sprintf(num2str(clock) + ".mat"));
    % load ('2025             12             10             12             11           0.16.mat'); 
    % load ('test_gen_oversampled_50_MHz.mat'); 

    figure(2);
    subplot(4,1,1)
    plot(s_to_subadc(1:100));
    subplot(4,1,2)
    plot(s_after_subadc);
    subplot(4,1,3)
    plot(adc_input_int(1:100,1));
    subplot(4,1,4)
    snr(s, sim_options.Fs, 2, 'aliased');
end