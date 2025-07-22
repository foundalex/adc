% 1) Джиган В.И, Адаптивные фильтры
% 2) Айфичер Э, Джервис Б, Цифровая обработка сигналов. Практический подход
% 3) Hu.M, Yi.P, (2022), Digital Calibration for Gain, Time Skew, and Bandwidth Mismatch 
% in Under-Sampling Time-Interleaved System
% 4) Behrouz Farhang-Boroujeny, Adaptive Filters Theory and Applications 

clear all;
clc;

%% Parameters;
N = 73;
Fs = 2000000000;           % samples per second
dt = 1/Fs;                  % seconds per sample
StopTime = 0.0001;         % seconds
t = (0:dt:StopTime-dt)';    % seconds
M = 2;                      % Num of sub ADC
%% Sine wave:
freq0 = 0;
for num = 5:5
    count = 5000000;
    freq0 = freq0 + count;  % frequency of fundamental tone
    freq1 = 2000000;

    % Main signal
    s0 = cos(2*pi*freq0*t);
    s1 = cos(2*pi*freq1*t);

    s = s0;% + s1;
    noise = awgn(s, 62);

    s = s + noise;

% ADC
% for k = 1:4000
%     adc_input(k,1) = s(k*M+1); % ADC0
%     adc_input(k,2) = s(k*M+2); % ADC1
%     adc_input(k,3) = s(k*M+3); % ADC2
%     adc_input(k,4) = s(k*M+4); % ADC3
% end

% ADC0
% s(1) - 1 Fs
% s(2) - 0.1 Fs
% s(3) - 0.2 Fs
% s(4) - 0.3 Fs
% s(5) - 0.4 Fs
% s(6) - 0.5 Fs
% s(7) - 0.6 Fs
% s(8) - 0.7 Fs
% s(9) - 0.8 Fs
% s(10) - 0.9 Fs
% s(11) - 1 Fs
% ADC1
% s(12) - 1 Fs

adc_input(:,1) = s(1:M*11:end-2*M*11);      % ADC0
adc_input(:,2) = s(1*11+1:M*11:end-2*M*11); % ADC1
% adc_input(:,3) = s(2*11+1:M*11:end-M*11); % ADC2
% adc_input(:,4) = s(3*11+1:M*11:end-M*11); % ADC3
% 
% adc_input(:,5) = s(4*11+1:M*11:end-M*11); % ADC4
% adc_input(:,6) = s(5*11+1:M*11:end-M*11); % ADC5
% adc_input(:,7) = s(6*11+1:M*11:end); % ADC6
% adc_input(:,8) = s(7*11+1:M*11:end); % ADC7


x_to_subadc = zeros(M*length(adc_input(:,1)),1);
x_to_subadc(1:M:end) = adc_input(:,1); 
x_to_subadc(2:M:end) = adc_input(:,2);
% x_to_subadc(3:M:end) = adc_input(:,3); 
% x_to_subadc(4:M:end) = adc_input(:,4); 
% 
% x_to_subadc(5:M:end) = adc_input(:,5); 
% x_to_subadc(6:M:end) = adc_input(:,6);
% x_to_subadc(7:M:end) = adc_input(:,7); 
% x_to_subadc(8:M:end) = adc_input(:,8); 

% adc_input(:,1) = adc_input(:,1); % ADC0
% adc_input(:,2) = adc_input(:,2) + noise1; % ADC1
% adc_input(:,3) = adc_input(:,3) + noise2; % ADC2
% adc_input(:,4) = adc_input(:,4) + noise3; % ADC3

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



% create main signal from adc
x_after_subadc = zeros(M*length(adc_input(:,1)),1);
x_after_subadc(1:M:end) = adc_input(:,1); 
x_after_subadc(2:M:end) = adc_input(:,2);
% x_after_subadc(3:M:end) = adc_input(:,3); 
% x_after_subadc(4:M:end) = adc_input(:,4); 
% 
% x_after_subadc(5:M:end) = adc_input(:,5); 
% x_after_subadc(6:M:end) = adc_input(:,6);
% x_after_subadc(7:M:end) = adc_input(:,7); 
% x_after_subadc(8:M:end) = adc_input(:,8); 

% plot(x_after_subadc(1:250))
% plot([adc_input(2:51,1), adc_input(1:50,2)]);

    % spectrumScope = spectrumAnalyzer(SampleRate=Fs/11, ...            
    %             AveragingMethod='exponential',ForgettingFactor=0.99, ...
    %             YLimits=[-30 10],ShowLegend=true, FrequencySpan="start-and-stop-frequencies", StartFrequency=0, StopFrequency=Fs/11/2);
    % spectrumScope(x_after_subadc);

    % figure(1);
    % subplot(3,1,1);
    % sfdr(x_to_subadc);
    % % title('SNR')
    % % xlabel('Нормированная частота') 
    % % ylabel('SNR (db)') 
    % % legend('до калибровки','после калибровки')
    % subplot(3,1,2);
    % sfdr(x_after_subadc);

    % figure(2);
    % subplot(3,1,1);
    % snr(x_to_subadc);
    % % title('SNR')
    % % xlabel('Нормированная частота') 
    % % ylabel('SNR (db)') 
    % % legend('до калибровки','после калибровки')
    % subplot(3,1,2);
    % snr(x_after_subadc);
    % % title('SFDR (db)')
    % % xlabel('Нормированная частота') 
    % % ylabel('SFDR (db)') 
    % % legend('до калибровки','после калибровки')
%% 
% Hu.M, Yi.P, (2022), Digital Calibration for Gain, Time Skew, and Bandwidth Mismatch 
% in Under-Sampling Time-Interleaved System

% fractional delay
% ADC0 - 1;
% ADC1 - 0.75
% ADC2 - 0.5
% ADC3 - 0.25
% for i = 1:M
%     delay_adc(i) = i/M; % (стр.6,(16))
% end
% delay_adc = fliplr(delay_adc);

n = (1:73);
ww = blackman(N);
% wvtool(blackman(N));


% adc_input_id(:,1) = adc_input(:,1);
% создаем массив [1:3] из задержанного сигнала ADC0 на различные значения [0.75 0.5 0.25]

delay_adc = [0.5; 0.5]; 

for i = 1:M
    % Method 1
    % % % delay_adc(i) = i/M; % (стр.6,(16))
    % h = sinc(n - (N - 1) / 2 - delay_adc(i));
    % te = h.' .* ww;
    % te = te ./ sum(te); 
    % freqz(te,1);
    % yri(:,i) = filter(te, 1, adc_input(:,1)); % (стр.6 (15))

    % Method 2
    [hri(:,i), i0(i), bw(i)] = designFracDelayFIR(delay_adc(i), N); 
    % plot([h.', hri(:,i)])
    [H1,w] = freqz(hri(:,i),1);
    plot(w/pi,mag2db(abs([H1])))
    % total_delay_fir(i) = i0(i) + delay_adc(i);
    fdfir = dsp.FIRFilter(hri(:,i).'); 

    yri(:,i) = fdfir(adc_input(:,i)); % (стр.6 (15))
    % plot(yri(:,i))

    plot([adc_input(1:150,2), yri(1:150,i)]);
    zz = yri((N-1)/2:end,i);
    adc_input_id(:,i) = zz;
    plot([adc_input(1:100,2), zz(1:100)]);
 
end
plot([adc_input(1:250,2), adc_input_id(1:250,2)]);

adc_input_id(9053:end,:) = [];
% % % aa = adc_input(1:end-((N-1)/2)+2,1);

% aa = [adc_input_id((N-1)/2:end-1,1)];
aa = adc_input_id(1:end,1);
bb = [adc_input_id(2:end,2); 0];
% cc = adc_input_id(1:end,3);
% dd = adc_input_id(1:end,4);
% dd = adc_input_id(1:end,4);
% figure(2);
% plot([adc_input(1:100,1), aa(1:100)]);
plot([adc_input(1:100,1), bb(1:100)]); %, cc(1:100), dd(1:100)]);


tt = 2;
x_after_subadc1 = zeros(tt*length(adc_input_id(:,1)),1);
% % x_after_subadc1(4:tt:end) = [0; adc_input(1:length(adc_input_id(1:end-1,1)),1)]; 
% 
ab = adc_input_id(:,2);
x_after_subadc1(2:tt:end) = ab; 
%
ab1 = [0; adc_input(1:end-((N-1)/2)-2,2)];
% x_after_subadc1(2:tt:end) = ab1; 
% 
x_after_subadc1(1:tt:end) = [adc_input_id(:,1)];
% % x_after_subadc1(3:tt:end) = adc_input_id(:,2);
% % x_after_subadc1(4:tt:end) = [0; adc_input_id(1:end-1,1)];

plot([ab1(1:100), ab(1:100)]);

    % spectrumScope = spectrumAnalyzer(SampleRate=Fs/11/M, ...            
    %             AveragingMethod='exponential',ForgettingFactor=0.99, ...
    %             YLimits=[-30 10],ShowLegend=true, FrequencySpan="start-and-stop-frequencies", StartFrequency=0, StopFrequency=Fs/11/4);
    % spectrumScope(ab1,ab);

% plot([adc_input_id(1:100,2), ab(1:100)])

% figure(3);


a1 = adc_input(1:9052,1);
a2 = adc_input(1:9052,2);
a3 = adc_input_id(:,1);
a4 = adc_input_id(:,2);
a5 = x_after_subadc(1:18104);
a6 = x_after_subadc1(1:end);
% plot([a1, a2]);


    figure(2);
    subplot(4,1,1);
    sfdr(a1, Fs/11/2);
    subplot(4,1,2);
    sfdr(a2, Fs/11/2);
    subplot(4,1,3);
    sfdr(a3, Fs/11/2);
    subplot(4,1,4);
    sfdr(a4, Fs/11/2);

    figure(3);
    subplot(2,1,1);
    sfdr(a5, Fs/11);
    subplot(2,1,2);
    sfdr(a6, Fs/11);


    % sfdr(a1, Fs/11);
    % subplot(4,1,4);
    % sfdr(a2, Fs/11);
    % 
    % figure(1);
    % subplot(4,1,1);
    % sfdr(adc_input(:,1));
    % subplot(4,1,2);
    % sfdr(adc_input(:,2));
    % % subplot(8,1,3);
    % % sfdr(adc_input(:,3));
    % % subplot(8,1,4);
    % % sfdr(adc_input(:,4));
    % subplot(4,1,3);
    % sfdr(adc_input_id(:,1));
    % subplot(4,1,4);
    % sfdr(adc_input_id(:,2));
    % % subplot(8,1,7);
    % % sfdr(adc_input_id(:,3));
    % % subplot(8,1,8);
    % % sfdr(adc_input_id(:,4));

    % spectrumScope = spectrumAnalyzer(SampleRate=Fs/11, ...            
    %             AveragingMethod='exponential',ForgettingFactor=0.99, ...
    %             YLimits=[-30 10],ShowLegend=true, FrequencySpan="start-and-stop-frequencies", StartFrequency=0, StopFrequency=Fs/11/4);
    % spectrumScope([x_after_subadc(1:8192),x_after_subadc1(1:8192)]);

% adc_input(:,1) = [adc_input(N/2-1:end,1); zeros(round(N/2)-2,1)];
% plot([adc_input(N/2:20+N/2-1,1), adc_input_id(1:20,1), adc_input_id(1:20,2), adc_input_id(1:20,3)]);
% plot([adc_input(3:21,1), adc_input_id(2:20,1), adc_input_id(2:20,2), adc_input_id(2:20,3)]);

% проверяем правильно ли мы задержали компоненты сигнала
% должен получиться синус как на входе фильтра
% aa = adc_input(N/2:end-(N/2),1); 
% aa = adc_input(3:end-(N/2)+2,1); 
% x_after_subadc1 = zeros(4*length(aa),1);
% x_after_subadc1(1:4:end) = aa; 
% x_after_subadc1(1:4:end) = adc_input_id(2:end,1);
% x_after_subadc1(2:4:end) = adc_input_id(2:end,2);
% x_after_subadc1(4:4:end) = adc_input_id(2:end,3);
% plot([x_after_subadc(1:250), x_after_subadc1(1:250)])

%% Least Mean algorithm

    % filter and coeff estimate first N words
    y_out1 = [];
    x3 = [];

    % считаем сигналы для ADC1-ADC3, т.к для ADC0 сигнал известен
    for z = 2:M
        y_out = 0;

        % создаем матрицу входного сигнала
        for i = 1:N
            x3(i,:) = adc_input(i:N+i-1,z).'; % (стр.6, (20))
        end

        % рассчитываем первые N коэффициентов адаптивного фильтра
        % сравнивая с задержанным сигналом ADC0 (adc_input_id)
        w1 = (x3'*x3) \ x3' * adc_input_id(1:N,z); % (стр.6, (19))
        % dd = det(x3);
        % ee = cond(x3);
        % w1 = inv(x3'*x3) * x3' * s(1:N);
        % w1 = lsqr(x3, adc_input_id(1:N,z));

        % умножаем входные слова на рассчитанные коэффициенты
        for k = 1:N
            y_out = y_out + w1(k)*adc_input(k,z); % (стр 5, (13))
            % y_out = w1(1)*x(j+1) + w1(2)*x(j+2) + w1(3)*x(j+3) +
            % w1(4)*x(j+4); % Behrouz Farhang-Boroujeny, Adaptive Filters
            % Theory and Applications  (стр. 414)
        end
        y_out1(1) = y_out;

    %% other wors
        for j = 1:length(adc_input(:,1))-2*N
    
            % shift to left matrix input signal. Refresh matrix input signal
            % for every new word
            for i = 1:N
                x3(i,:) = [x3(i,2:N), 0];
                x3(i,N) = adc_input(j+N-1+i,z); % (стр.6, (20))
            end

            % estimate coeff
            w1 = ((x3'*x3) \ x3') * adc_input_id(j+1:N+j,z); % (стр.6, (19))
            % w1 = lsqr(x3, adc_input_id(j+1:N+j,z));

            % filter input signal. Mult input words on coeff
            y_out = 0;
            for k = 1:N
                y_out = y_out + w1(k)*adc_input(j+k,z); % (стр 5, (13))
            end
            y_out1(j+1) = y_out;
        end


        y_out11(1:j+1,z-1) =  y_out1(1:j+1).';
        error_out(1:j+1,z-1) = y_out1(1:j+1).' ./ adc_input_id(1:j+1,z);

        figure(5);
        plot([adc_input(1:j+1,z), y_out1(1:j+1).']);
        % plot(error_out);

    end

    % create main signal from sub-adc
    x_after_adc = zeros(j*M+M,1);
    x_after_adc(4:M:end) = adc_input_id(1:j+1,1); 
    x_after_adc(3:M:end) = y_out11(1:j+1,1); 
    x_after_adc(2:M:end) = y_out11(1:j+1,2); 
    x_after_adc(1:M:end) = y_out11(1:j+1,3); 

    % x_after_adc(5:M:end) = y_out11(1:j+1,4); 
    % x_after_adc(6:M:end) = y_out11(1:j+1,5); 
    % x_after_adc(7:M:end) = y_out11(1:j+1,6); 
    % x_after_adc(8:M:end) = y_out11(1:j+1,7); 

    % plot(x_after_adc(1:50));

    figure(3);
    subplot(2,1,1)
    plot([x_after_subadc(1:500), x_after_adc(1:500)])
    % title('Отношение между отсчетами I-составляющей')
    xlabel('Номер отсчета') 
    ylabel('Амплитуда') 
    legend({'Исходный сигнал','Сигнал с выхода адаптивного фильтра'},'Location','northeast')
    subplot(2,1,2)
    plot([error_out(:,1)]); %, error_out(:,2), error_out(:,3)]);
    title('Относительная ошибка между исходным сигналом и выходом адаптивного фильтра')
    xlabel('Номер отсчета') 
    ylabel('Отношение') 

    % spectrumScope = spectrumAnalyzer(SampleRate=Fs/11, ...            
    %             AveragingMethod='exponential',ForgettingFactor=0.99, ...
    %             YLimits=[-30 10],ShowLegend=true, FrequencySpan="start-and-stop-frequencies", StartFrequency=0, StopFrequency=Fs/11/2);
    % spectrumScope([x_after_subadc(1:M*nn),x_after_adc]);
%% Measurements
    snr_input(num+1) = snr(x_after_subadc(1:17596));
    snr_output(num+1) = snr(x_after_adc(1:end));
    sfdr_input(num+1) = sfdr(x_after_subadc(1:17596), Fs/11/M);
    sfdr_output(num+1) = sfdr(x_after_adc(1:end), Fs/11/M);
    norm_freq(num+1) = freq0/(Fs/11/M);

    figure(4);
    subplot(2,1,1)
    plot(norm_freq, snr_input, '-o', norm_freq, snr_output, '-o');
    title('SNR')
    xlabel('Нормированная частота') 
    ylabel('SNR (db)') 
    legend('до калибровки','после калибровки')
    subplot(2,1,2)
    plot(norm_freq, sfdr_input, '-o', norm_freq, sfdr_output, '-o');
    title('SFDR (db)')
    xlabel('Нормированная частота') 
    ylabel('SFDR (db)') 
    legend('до калибровки','после калибровки')

    % figure(1);
    % subplot(3,1,3);
    % sfdr(x_after_adc);
    % 
    % figure(2);
    % subplot(3,1,3);
    % snr(x_after_adc);

end


