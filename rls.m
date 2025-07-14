% 1) Джиган В.И, Адаптивные фильтры
% 2) Айфичер Э, Джервис Б, Цифровая обработка сигналов. Практический подход
% 3) Hu.M, Yi.P, (2022), Digital Calibration for Gain, Time Skew, and Bandwidth Mismatch 
% in Under-Sampling Time-Interleaved System
% 4) Behrouz Farhang-Boroujeny, Adaptive Filters Theory and Applications 

clear all;
clc;
%% Time specifications:
freq = 1000000;             % frequency of fundamental tone
N = 73;
Fs = 20000000;              % samples per second
dt = 1/Fs;                  % seconds per sample
StopTime = 0.001;           % seconds
t = (0:dt:StopTime-dt)';    % seconds
%% Sine wave:

% Main signal
s = 0.07*cos(2*pi*freq*t);
noise = awgn(s, 60);
s = s + noise;

% Num of sub ADC
M = 4;

% ADC
for k = 1:4000
    adc_input(k,1) = s(k*M+1); % ADC0
    adc_input(k,2) = s(k*M+2); % ADC1
    adc_input(k,3) = s(k*M+3); % ADC2
    adc_input(k,4) = s(k*M+4); % ADC3
end

adc_input(4000+1:end,:) = [];


% noise1 = awgn(adc_input(:,1), 61);
% adc_input(:,1) = adc_input(:,1)+noise1; % ADC0
% adc_input(:,2) = adc_input(:,2)+noise1; % ADC1
% adc_input(:,3) = adc_input(:,3)+noise1; % ADC2
% adc_input(:,4) = adc_input(:,4)+noise1; % ADC3


% % model offset error
% adc_input(:,2) = adc_input(:,2) + 0.002;
% adc_input(:,3) = adc_input(:,3) + 0.05;
% adc_input(:,4) = adc_input(:,4) + 0.010;
% 
% model gain error
% adc_input(:,2) = adc_input(:,2) * 1.4;
% adc_input(:,3) = adc_input(:,3) * 1.2;
% adc_input(:,4) = adc_input(:,4) * 1.5;

% model timing skew
adc_input(:,2) = adc_input(:,3); % 0.5/fs

% create main signal from adc
x_after_subadc = zeros(M*length(adc_input(:,1)),1);

x_after_subadc(1:M:end) = adc_input(:,1); 
x_after_subadc(2:M:end) = adc_input(:,2);
x_after_subadc(3:M:end) = adc_input(:,3); 
x_after_subadc(4:M:end) = adc_input(:,4); 
% 
plot(x_after_subadc(1:250))
%% 
% Hu.M, Yi.P, (2022), Digital Calibration for Gain, Time Skew, and Bandwidth Mismatch 
% in Under-Sampling Time-Interleaved System
for i = 1:M-1
    delay_adc(i) = i/M; % (стр.6,(16))
   
end
delay_adc = fliplr(delay_adc);

for i = 1:M-1
    hri(:,i) = designFracDelayFIR(delay_adc(i),N); 
    yri(:,i) = filter(hri(:,i),1,adc_input(:,1)); % (стр.6 (15))
    % plot([adc_input(1:200,1), yri(N/2:200+(N/2)-1)]);
    adc_input(1:end-N/2+1,i+1) = yri(N/2:end,i); % (N/2 (delete transient processes), стр. 6, (21))
    % plot([adc_input(1:200,1), adc_input(1:200,2)]);
end

adc_input(end-N/2+1:end,:) = [];
plot([adc_input(1:20,1), adc_input(1:20,2), adc_input(1:20,3), adc_input(1:20,4)]);

% create main signal from adc
x_after_subadc1 = zeros(4*length(adc_input(:,1)),1);

x_after_subadc1(1:M:end) = adc_input(:,1); 
x_after_subadc1(2:M:end) = adc_input(:,2);
x_after_subadc1(3:M:end) = adc_input(:,3);
x_after_subadc1(4:M:end) = adc_input(:,4);
plot([x_after_subadc(1:250), x_after_subadc1(1:250)])

%% Least Mean algorithm

    % filter and coeff estimate first N words
    y_out1 = [];
    x3 = [];

    for z = 2:M
        y_out = 0;

        % create matrix for input signal
        for i = 1:N
            x3(i,:) = adc_input(i:N+i-1,z).'; % (стр.6, (20))
        end

        % calculate first N coefficients
        w1 = (x3'*x3) \ x3' * adc_input(1:N,1); % (стр.6, (19))
        % dd = det(x3);
        % ee = cond(x3);
        % w1 = inv(x3'*x3) * x3' * s(1:N);
        % st = lsqr(x3, s(1:N));

        % filter first N words
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
            % w1 = inv(x3'*x3) * x3' * x_adc0(j+1:N+j);
            w1 = ((x3'*x3) \ x3') * adc_input(j+1:N+j,z); % (стр.6, (19))
            % st = lsqr(x3, x_adc0(j+1:N+j));

            % filter input signal. Mult input words on coeff
            y_out = 0;
            for k = 1:N
                y_out = y_out + w1(k)*adc_input(j+k,z); % (стр 5, (13))
            end
            y_out1(j+1) = y_out;
        end

        nn = 3819;
        y_out11(1:nn,z-1) =  y_out1(1:nn).';
        error_out(1:nn,z-1) = y_out1(1:nn).' ./ adc_input(1:nn,z);
    end

    % create main signal from sub-adc
    x_after_adc = zeros(nn*M,1);
    x_after_adc(1:4:end) = adc_input(1:nn,1); 
    x_after_adc(2:4:end) = y_out11(1:nn,1); 
    x_after_adc(3:4:end) = y_out11(1:nn,2); 
    x_after_adc(4:4:end) = y_out11(1:nn,3); 

    subplot(2,1,1)
    plot([x_after_subadc(1:nn*M), x_after_adc])
    % title('Отношение между отсчетами I-составляющей')
    xlabel('Номер отсчета') 
    ylabel('Амплитуда') 
    legend({'Исходный сигнал','Сигнал + Шум','Сигнал с выхода адаптивного фильтра'},'Location','northeast')
    subplot(2,1,2)
    plot([error_out(:,1)]);% error_out(:,2), error_out(:,3)])
    title('Относительная ошибка между исходным сигналом и выходом адаптивного фильтра')
    xlabel('Номер отсчета') 
    ylabel('Отношение') 

    spectrumScope = spectrumAnalyzer(SampleRate=Fs, ...            
                AveragingMethod='exponential',ForgettingFactor=0.99, ...
                YLimits=[-30 10],ShowLegend=true, FrequencySpan="start-and-stop-frequencies", StartFrequency=0, StopFrequency=Fs/2);
    spectrumScope([x_after_subadc(1:M*nn),x_after_adc]);