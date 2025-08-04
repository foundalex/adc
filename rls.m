% 1) Джиган В.И, Адаптивные фильтры
% 2) Айфичер Э, Джервис Б, Цифровая обработка сигналов. Практический подход
% 3) Hu.M, Yi.P, (2022), Digital Calibration for Gain, Time Skew, and Bandwidth Mismatch 
% in Under-Sampling Time-Interleaved System
% 4) Behrouz Farhang-Boroujeny, Adaptive Filters Theory and Applications 

close all;
clear all;
clc;

%% Parameters;
N = 73;
Fs = 2000000000;           % Fs = 2 GHz
dt = 1/Fs;                 % seconds per sample
StopTime = 0.0001;         % seconds
t = (0:dt:StopTime-dt)';   % seconds
M = 2;                     % Num of sub ADC
Fc = Fs/M;
%% Sine wave:
freq0 = 470000000; % 44 MHz
for num = 0:24
    count = 10000000;
    freq0 = freq0 + count;       % frequency of fundamental tone
    Z = ceil(freq0/(Fs/2/M));    % Nyquist zone
    % Main signal
    s0 = cos(2*pi*freq0*t);
    s = s0;% + s1;
    noise = awgn(s,65);
    s = s + noise;

    % ADC0
    % s(1) - 1 Fs, s(2) - 0.1 Fs, s(3) - 0.2 Fs, s(4) - 0.3 Fs
    % s(5) - 0.4 Fs, s(6) - 0.5 Fs, s(7) - 0.6 Fs, s(8) - 0.7 Fs
    % s(9) - 0.8 Fs, s(10) - 0.9 Fs, s(11) - 1 Fs
    % ADC1
    % s(12) - 1 Fs

    adc_input = [];
    for i = 1:M
        adc_input(:,i) = s(i:M:end);
    end 

	% adc_input(:,1) = s(1:M*11:end-2*M*11);        % ADC0
	% adc_input(:,2) = s(1*11+1:M*11:end-2*M*11);   % ADC1
	% adc_input(:,3) = s(2*11+1:M*11:end-M*11);     % ADC2
	% adc_input(:,4) = s(3*11+1:M*11:end-M*11);     % ADC3
	% 
	% adc_input(:,5) = s(4*11+1:M*11:end-M*11);     % ADC4
	% adc_input(:,6) = s(5*11+1:M*11:end-M*11);     % ADC5
	% adc_input(:,7) = s(6*11+1:M*11:end);          % ADC6
	% adc_input(:,8) = s(7*11+1:M*11:end);          % ADC7

	% add noise
	for i = 1:M
		noise1(:,i) = awgn(adc_input(:,i),60);
		adc_input(:,i) = adc_input(:,i) + noise1(:,i);
	end

	% исходный сигнал до искажений
	x_to_subadc = zeros(M*length(adc_input(:,1)),1);
	for i = 1:M
		x_to_subadc(i:M:end) = adc_input(:,i); 
	end

	% model timing skew
	% adc_input(:,2) = s(17:M*11:end-M*11); % 0.5/fs
	% adc_input(:,3) = s(25:M*11:end-M*11); % 0.2/fs
	% adc_input(:,4) = s(35:M*11:end-M*11); % 0.1/fs
	% adc_input(:,5) = s(49:M*11:end-M*11); % 0.4/fs
	% adc_input(:,6) = s(59:M*11:end-M*11); % 0.3/fs
	% adc_input(:,7) = s(68:M*11:end); % 0.1/fs
	% adc_input(:,8) = s(80:M*11:end); % 0.2/fs

	% % model offset error
	% adc_input(:,2) = adc_input(:,2) + 0.002;
	% adc_input(:,3) = adc_input(:,3) + 0.05;
	% adc_input(:,4) = adc_input(:,4) + 0.010;
	% 
	% model gain error
	% adc_input(:,2) = adc_input(:,2) * 1.4;
	% adc_input(:,3) = adc_input(:,3) * 1.3;
	% adc_input(:,4) = adc_input(:,4) * 1.2;
	% adc_input(:,5) = adc_input(:,5) * 1.1;
	% adc_input(:,6) = adc_input(:,6) * 1.2;
	% adc_input(:,7) = adc_input(:,7) * 1.3;
	% adc_input(:,8) = adc_input(:,8) * 1.1;

	% Мейн сигнал с искажениями
	x_after_subadc = zeros(M*length(adc_input(:,1)),1);
	for i = 1:M
		x_after_subadc(i:M:end) = adc_input(:,i); 
	end

    % figure(5)
    % plot([adc_input(1:100,1), adc_input(1:100,3)]);




	%% 
	% Hu.M, Yi.P, (2022), Digital Calibration for Gain, Time Skew, and Bandwidth Mismatch 
	% in Under-Sampling Time-Interleaved System

	% Дробная задержка сигнала ADC0
	[yri_cut] = fractional_delays(adc_input, M, N, Z);

	% yri_cut(:,1) = [zeros(36,1); yri_cut(1:end-36,1)];
	figure(10);
	plot([yri_cut(1:50,1), yri_cut(1:50,2)]); %, yri_cut(1:50,3), yri_cut(1:50,4)]);

	% сигнал после фильтров с задержками
	sig_adc = zeros(M*length(yri_cut(:,1)),1);
	for i = 1:M
		sig_adc(i:M:end) = yri_cut(:,M-i+1);
	end

	figure(11);
	plot([x_after_subadc(1:100), sig_adc(1:100)]); %((N-1)/2)*M
    % 
	figure(1);
	subplot(3,1,1);
	sfdr(x_to_subadc(1:length(sig_adc)), Fs);
	subplot(3,1,2);
	sfdr(x_after_subadc(1:length(sig_adc)), Fs);
	subplot(3,1,3);
	sfdr(sig_adc(1:length(sig_adc)), Fs);
	% 
	figure(2);
	subplot(3,1,1);
	snr(x_to_subadc(1:length(sig_adc)), Fs);
	subplot(3,1,2);
	snr(x_after_subadc(1:length(sig_adc)), Fs);
	subplot(3,1,3);
	snr(sig_adc(1:length(sig_adc)), Fs);

	spectrumScope = spectrumAnalyzer(SampleRate=Fs, ...            
	            AveragingMethod='exponential',ForgettingFactor=0, ...
	            YLimits=[-30 10],ShowLegend=true, Method='Welch');

	spectrumScope.WindowLength = 2048;
	spectrumScope.FrequencyResolutionMethod = "window-length";
	spectrumScope.PlotAsTwoSidedSpectrum=false;
	spectrumScope.DistortionMeasurements.Enabled = true;


	spectrumScope([sig_adc(4096:19000)]);
% end


    %% Least Mean algorithm
    % считаем сигналы для ADC1-ADC3, т.к для ADC0 сигнал известен
    for z = 2:M
        %% блок для расчета первых N коэффициентов фильтра
        y_out(1) = 0;
        % создаем матрицу входного сигнала
        for i = 1:N
            x3(i,:) = adc_input(i:N+i-1,z).'; % (стр.6, (20))
        end
        % рассчитываем первые N коэффициентов адаптивного фильтра
        % сравнивая с задержанным сигналом ADC0 (yri_cut)
        w1 = (x3'*x3) \ x3' * yri_cut(1:N,z); % (стр.6, (19))
        % w1 = lsqr(x3, adc_input_id(1:N,z));

        % умножаем входные слова на рассчитанные коэффициенты
        for k = 1:N
            y_out(1) = y_out(1) + w1(k)*adc_input(k,z); % (стр 5, (13))
            % y_out = w1(1)*x(j+1) + w1(2)*x(j+2) + w1(3)*x(j+3) +
            % w1(4)*x(j+4); % Behrouz Farhang-Boroujeny, Adaptive Filters
            % Theory and Applications  (стр. 414)
        end

        %% пересчет коэффициентов с приходом каждого слова
        for j = 1:length(yri_cut(:,1))-2*N
  
            % shift to left matrix input signal. Refresh matrix input signal
            % for every new word
            for i = 1:N
                x3(i,:) = [x3(i,2:N), 0];
                x3(i,N) = adc_input(j+N-1+i,z); % (стр.6, (20))
            end

            % estimate coeff
            w1 = ((x3'*x3) \ x3') * yri_cut(j+1:N+j,z); % (стр.6, (19))
            % w1 = lsqr(x3, adc_input_id(j+1:N+j,z));

            % filter input signal. Mult input words on coeff
            y_out = 0;
            for k = 1:N
                y_out = y_out + w1(k)*adc_input(j+k,z); % (стр 5, (13))
            end
            y_out1(j+1) = y_out;
        end

        y_out_array(1:j+1,z-1) =  y_out1(1:j+1).';
        error_out(1:j+1,z-1) = y_out1(1:j+1).' ./ yri_cut(1:j+1,z);

    end

    % create main signal after LS algorithm
    x_after_adc = zeros(length(y_out_array)*M,1);
    for i = 1:M
        if i == 1
            x_after_adc(i:M:end) = yri_cut(1:length(y_out1));
        else
            x_after_adc(i:M:end) = y_out_array(:,i-1);
        end
    end

    %% Measurements
    figure(4);
    subplot(2,1,1)
    plot([x_after_subadc(1:150), x_after_adc(1:150)])
    % title('Отношение между отсчетами I-составляющей')
    xlabel('Номер отсчета') 
    ylabel('Амплитуда') 
    legend({'Исходный сигнал','Сигнал с выхода адаптивного фильтра'},'Location','northeast')
    subplot(2,1,2)
    plot([error_out(:,1)]); %, error_out(:,2), error_out(:,3)]);
    title('Относительная ошибка между исходным сигналом и выходом адаптивного фильтра')
    xlabel('Номер отсчета') 
    ylabel('Отношение') 

    figure(1);
    subplot(4,1,1);
    sfdr(x_to_subadc(1:length(x_after_adc)), Fs);
    subplot(4,1,2);
    sfdr(x_after_subadc(1:length(x_after_adc)), Fs);
    subplot(4,1,3);
    sfdr(sig_adc(1:length(x_after_adc)), Fs);
    subplot(4,1,4);
    sfdr(x_after_adc(1:length(x_after_adc)), Fs);

    figure(2);
    subplot(4,1,1);
    snr(x_to_subadc(1:length(x_after_adc)), Fs);
    subplot(4,1,2);
    snr(x_after_subadc(1:length(x_after_adc)), Fs);
    subplot(4,1,3);
    snr(sig_adc(1:length(x_after_adc)), Fs);
    subplot(4,1,4);
    snr(x_after_adc(1:length(x_after_adc)), Fs);

    snr_input(num+1) = snr(x_after_subadc);
    snr_output(num+1) = snr(x_after_adc);
    sfdr_input(num+1) = sfdr(x_after_subadc, Fs);
    sfdr_output(num+1) = sfdr(x_after_adc, Fs);
    norm_freq(num+1) = freq0/(Fs/2);

    ss(num+1) = sfdr(sig_adc(1:length(sig_adc)));

end


figure(5);
subplot(3,1,1)
plot(norm_freq, snr_input, '-o', norm_freq, snr_output, '-o');
title('SNR')
xlabel('Нормированная частота') 
ylabel('SNR (db)') 
legend('до калибровки','после калибровки')
subplot(3,1,2)
plot(norm_freq, sfdr_input, '-o', norm_freq, sfdr_output, '-o');
title('SFDR (db)')
xlabel('Нормированная частота') 
ylabel('SFDR (db)') 
legend('до калибровки','после калибровки')
subplot(3,1,3)
plot(norm_freq, ss, '-o');