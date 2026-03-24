function [s, s_to_subadc, adc_input_double, adc_input_int, s_after_subadc, Z] = gen_oversampled_signal(sim_options)

    dt = 1/sim_options.Fs;                                                              % seconds per sample
    t = 0:dt:sim_options.StopTime;                                                      % seconds

    % Create main signal with noise in double


    s1 = 1*sin(2*pi*sim_options.freq*t);
    
    % save (sprintf(num2str(clock) + ".mat"));
    % load ('2026              3              6             15             40         18.288.mat'); 

    % sim_options.MODEL_ERROR = true;
    % sim_options.time_skew_array = [0.01, 0.04, 0.03, 0];

    % ss = -sin(2*pi*sim_options.freq*t);


    % s2 = 0.25*sin(2*pi*(1440000000)*t);
    % s3 = 0.025*sin(2*pi*(sim_options.freq+800000000)*t);
    % s4 = 0.025*sin(2*pi*(sim_options.freq+1500000000)*t);
    % s5 = 0.025*sin(2*pi*(sim_options.freq+700000000)*t);
    % s6 = 0.025*sin(2*pi*(sim_options.freq+900000000)*t);
    % s7 = 0.025*sin(2*pi*(sim_options.freq+1000000000)*t);
    % s8 = 0.025*sin(2*pi*(sim_options.freq+1400000000)*t);
    % s9 = 0.025*sin(2*pi*(sim_options.freq+1500000000)*t);
    % s10 = 0.025*sin(2*pi*(sim_options.freq+1700000000)*t);


    % figure(2);
    % plot([s1(1:100)', ss(1:100)']);

    %% АМ-модуляция
    % sim_options.Fs = 1000000000; % частота дискретизации
    % dt = 1/sim_options.Fs;                                                              % seconds per sample
    % t = 0:dt:sim_options.StopTime; 
    % s1 = 1+0.5*cos(2*pi*sim_options.freq*t); % информационный сигнал (постоянная составляющая и коэффициент модуляции)
    % sm = cos(2*pi*200000000*t); % несущее колебание 
    % % частота несущего колебания должна быть >> информационного сигнала
    % en = s1.*sm;
    % 
    % figure(77);
    % subplot(2,1,1)
    % plot([en(1:600)]);
    % subplot(2,1,2)
    % snr(en, sim_options.Fs/1);
    %%

    % figure(3);
    % plot([s1(1:100)', s2(1:100)']);

    s = s1; %+ s2; % + s3; % + s4 + s5 + s6 + s7 + s8 + s9 + s10;

    % Sine1 = dsp.SineWave(Frequency=sim_options.freq,SampleRate=sim_options.Fs,SamplesPerFrame=2*10^6);
    % y = Sine1();
    % figure(4);
    % plot([y(1:100),s(1:100)'])


    s = awgn(s,sim_options.SNR(1));
    % y = awgn(y,sim_options.SNR(1));
    % s_after_subadc = y;
    % s = awgn(s, 60);
    % s = s + noise;
    % [pxx,f] = periodogram(s); 

    step = sim_options.M*sim_options.Inter;
    begin = (1:sim_options.M).*sim_options.Inter-(sim_options.Inter-1);

    if sim_options.MODEL_ERROR == true
        dlin = floor((length(s)-(sim_options.time_skew_array(1:sim_options.M-1).*sim_options.Inter))./(sim_options.M*sim_options.Inter));
        dlin = dlin - 1;
    else
        dlin = (floor((length(s)-begin)./(sim_options.M*sim_options.Inter)));
    end

    % adc_input = zeros(max(dlin),sim_options.M);
    % adc_input_int = int16(adc_input);

    ended = max(dlin)*sim_options.Inter*sim_options.M;
    
    %% разбиваем входной сигнал на сигналы для суб-АЦП
    
    for i = 1:sim_options.M  
        % если ошибки суб-АЦП включены, то добавляем time s
        adc_input_good(:,i) = s(begin(i):step:ended);
        if sim_options.MODEL_ERROR == true
            % time skew
            adc_input(:,i) = s(begin(i) + sim_options.time_skew_array(i)*sim_options.Inter:step:ended);
            % adc_input_sin(:,i) = ss(begin(i) + sim_options.time_skew_array(i)*sim_options.Inter:step:ended);
            % offset
            adc_input(:,i) = adc_input(:,i) + sim_options.offset_error_array(i);
            % gain
            adc_input(:,i) = adc_input(:,i) * sim_options.gain_error_array(i);
        else
            if i == 1
                adc_input(:,i) = s(begin(i):step:end-1);
            else
                adc_input(:,i) = s(begin(i):step:end);
            end
            % adc_input_sin(:,i) = ss(begin(i):step:ended);
        end
    end




    %% добавляем к каждому суб-АЦП шум
    adc_input_double = adc_input;
    % for i = 1:sim_options.M 
    %     adc_input(:,i) = awgn(adc_input(:,i), sim_options.SNR(i) , "measured");
    %     adc_input_double(:,i) = awgn(adc_input_double(:,i), sim_options.SNR(i) , "measured");
    %     adc_input_good(:,i) = awgn(adc_input_good(:,i), sim_options.SNR(i) , "measured");
    % end

    for i = 1:sim_options.M 
        s_fi(:,i) = fi(adc_input(:,i),1,12,11);
        s_int(:,i) = int16(round(s_fi(:,i)*2^11));
    end

    % перевод в инты
    for i = 1:sim_options.M
        adc_input_int(:,i) = int16(round(fi(adc_input(:,i),1,12,11)*sim_options.Bit));
        adc_input_int_good(:,i) = int16(round(fi(adc_input_good(:,i),1,12,11)*sim_options.Bit));
    end


    % e0 = zeros(1500,1);
    % e1 = ones(1500,1);
    % ee = [e0; e1];
    % 
    % spectrumScope = spectrumAnalyzer(SampleRate=sim_options.Fs, ...            
    %         AveragingMethod='exponential',ForgettingFactor=0.99, ...
    %         YLimits=[-30 10],ShowLegend=true);
    % 
    % spectrumScope(ee);
    % figure(33); plot(ee);


	% исходный сигнал до искажений
	s_to_subadc = zeros(sim_options.M*length(adc_input_int_good(:,1)),1);
	for i = 1:sim_options.M
		s_to_subadc(i:sim_options.M:end) = adc_input_int_good(:,i); 
    end

    % сигнал после искажений
    s_after_subadc = zeros(sim_options.M*length(adc_input(:,1)),1);
	for i = 1:sim_options.M
        s_after_subadc(i:sim_options.M:end) = adc_input_int(:,i); 
    end

    Z = ceil(sim_options.freq/(sim_options.Fs/sim_options.Inter/2/sim_options.M));      % Nyquist zone

    % save (sprintf(num2str(clock) + ".mat"));
    % load ('2026              2             20             12              0         14.172.mat'); 
    % load ('test_gen_oversampled_50_MHz.mat'); 

    % [pxx,f] = periodogram(s_after_subadc);

    % figure(2);
    % plot([adc_input(1:100,1), adc_input(1:100,2)]);

    figure(3);
    subplot(4,1,1)
    plot(s_to_subadc(1:100));
    subplot(4,1,2)
    plot(s_after_subadc(1:100));
    % subplot(5,1,3)
    % plot(adc_input_int(1:100,1));
    subplot(4,1,3)
    snr(s_to_subadc, sim_options.Fs/sim_options.Inter);
    subplot(4,1,4)
    snr(s_after_subadc, sim_options.Fs/sim_options.Inter);

    figure(4);
    for i = 1:sim_options.M
        subplot(sim_options.M,1,i)
        plot(adc_input_int(1:100,i));
        title(['Выход АЦП' num2str(i)])
        xlabel('Амплитуда') 
        ylabel('Номер отсчета') 
    end

end