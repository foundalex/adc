clear all;
clc;

% function [n,x,s,fs,s1] = genData(numPts, freq, filt, nVar, SNR)
%     % Generate time values
%     t = linspace(0,1,numPts)';
%     fs = numPts;
% 
%     % Generate tone
%     s = sin(2*pi*freq*t);
%     s1 = sin(2*pi*freq*(t+3.4e4));
%     % Generate noise
%     n = sqrt(nVar)*randn(numPts,1);
% 
%     % Filter noise
%     addnoise = filter(filt, 1, n);
% 
%     % Plot filter
%     freqz(filt,1,1000)
% 
%     % Adjust SNR of tone
%     s = s/sqrt(var(s)/(10^(SNR/10)*var(n)));
%     disp(['Calculated SNR = ' num2str(10*log10(var(s)/var(n)))])
% 
%     % Add noise to signal
%     x = s + addnoise;
% end

% % Generate the data!
% [n,x,s,fs,s1] = genData(numPts, freq, filt, nVar, SNR);
% % parameter estimator menggunakan least square
% 
% Y = ([1 : 100] + rand(1, 1))'; % proses markov orde 1
% X = ([100 : 199] + rand(1, 1))';
% st_time = tic();
% % jumlah parameter
% np = 2;
% % orde model
% nmdl = np / 2;
% phi = [];
% N = 40;
% for m = nmdl + 1 : N - nmdl
%     phi(m, :) = [-Y(m-1:-1:m-nmdl, :)' X(m-1:-1:m-nmdl, :)'];
% end
% Y2 = Y(nmdl + 1 : N);
% % Parameter yang diestimasi
% beta = inv(phi'*phi) * phi' * Y2;
% Y3 = phi * beta;
% err = Y2 - Y3;
% emse = mse(err);
% t = 1 : length(Y3);
% figure;
% plot(t, Y2, t, Y3);
% title('Data Asli vs Estimasi');
% figure;
% plot(t, err);
% title('Galat');
% cnsmd_time = toc(st_time)

freq = 1000;             % frequency of fundamental tone
N = 73;
%% Time specifications:
Fs = 20000;                 % samples per second
dt = 1/Fs;                   % seconds per sample
StopTime = 0.5;             % seconds
t = (0:dt:StopTime-dt)';     % seconds
%% Sine wave:
s = 0.07*cos(2*pi*freq*t);
noise = awgn(s, 60);
s = s + noise;
M = 2;

% ADC
for k = 1:2000
    s0(k,:) = s(k*M+1); % ADC0
    s1(k,:) = s(k*M+2); % ADC1
end

%offset
s1 = s1 + 0.002;

%gain
s1 = s1 * 1.4;

plot([s0,s1]);

x_adc0 = s0;% + noise;
x_adc1 = s1;% + noise;

% create main signal from adc
x_after_subadc = zeros(2*length(s0),1);
x_after_subadc(1:2:end) = x_adc0; 
x_after_subadc(2:2:end) = x_adc1; 
plot(x_after_subadc(1:250))

spectrumScope = spectrumAnalyzer(SampleRate = Fs, ...            
            AveragingMethod='exponential',ForgettingFactor=0.99, ...
            YLimits=[-30 10],ShowLegend=true);
% spectrumScope([x_adc0, x_adc1]);
spectrumScope([x_after_subadc]);


for z = 1:1

    y_out1 = [];
    x3 = [];
    y_out = 0;

    % create matrix for input signal
    for i = 1:N
        x3(i,:) = x_adc1(i:N+i-1).';
    end

    w1 = (x3'*x3) \ x3' * x_adc0(1:N);
    % dd = det(x3);
    % ee = cond(x3);
    % w1 = inv(x3'*x3) * x3' * s(1:N);
    % st = lsqr(x3, s(1:N));

    for k = 1:N
        % y_out = w1(1)*x(j+1) + w1(2)*x(j+2) + w1(3)*x(j+3) + w1(4)*x(j+4);
        y_out = y_out + w1(k)*x_adc1(k);
    end
    % y_out = w1(1)*x(1) + w1(2)*x(2) + w1(3)*x(3) + w1(4)*x(4); % + st(5)*x(5) + st(6)*x(6) + st(7)*x(7) + st(8)*x(8) + st(9)*x(9) + st(10)*x(10);
    y_out1(1) = y_out;

    % work
    % for i = 1:2
    %   x3 = [x(i), x(i+1); x(i+1), x(i+2)];
    %   y3 = [s(i); s(i+1)];
    %   w1 = inv(x3'*x3) * x3' * y3;
    %   st = lsqr(x3, y3);
    %   y = st(1)*x(i)+st(2)*x(i+1);
    % end
%%

    for j = 1:length(x_adc0)-2*N
    
        % shift to left matrix input signal
        for i = 1:N
            x3(i,:) = [x3(i,2:N), 0];
            x3(i,N) = x_adc1(j+N-1+i);
        end

        % estimate coeff
        % w1 = inv(x3'*x3) * x3' * x_adc0(j+1:N+j);
        w1 = ((x3'*x3) \ x3') * x_adc0(j+1:N+j);
        % st = lsqr(x3, x_adc0(j+1:N+j));

        y_out = 0;
        % filter input signal
        for k = 1:N
            % y_out = w1(1)*x(j+1) + w1(2)*x(j+2) + w1(3)*x(j+3) + w1(4)*x(j+4);
            y_out = y_out + w1(k)*x_adc1(j+k);
        end

        y_out1(j+1) = y_out;

    end

    nn = 1855;

    y_out11(1:nn,z) =  y_out1(1:nn).';
    error_out(1:nn,z) = y_out1(1:nn).' ./ x_adc0(1:nn);

end

    % subplot(2,1,1)
    % plot([x_adc0(1:nn), x_adc1(1:nn), y_out11(1:nn)])
    % % title('Отношение между отсчетами I-составляющей')
    % xlabel('Номер отсчета') 
    % ylabel('Амплитуда') 
    % legend({'Исходный сигнал','Сигнал + Шум','Сигнал с выхода адаптивного фильтра'},'Location','northeast')
    % subplot(2,1,2)
    % plot([error_out(:,1)]);% error_out(:,2), error_out(:,3)])
    % title('Относительная ошибка между исходным сигналом и выходом адаптивного фильтра')
    % xlabel('Номер отсчета') 
    % ylabel('Отношение') 

    % create main signal from sub-adc
    x_after_adc = zeros(nn*2,1);

    x_after_adc1 = zeros(nn,1);
    x_after_adc2 = zeros(nn,1);

    for ty = 1:nn
        a = 1+((ty - mod(ty,2))/2);
        x_after_adc1(ty) = x_adc0(1+(ty - mod(ty,2))/2);
        x_after_adc2(ty) = y_out11(1+(ty - mod(ty,2))/2);
    end

    x_after_adc(1:2:end) = x_adc0(1:nn); 
    x_after_adc(2:2:end) = y_out11; 
    plot(x_after_adc)

    spectrumScope = spectrumAnalyzer(SampleRate=Fs, ...            
                AveragingMethod='exponential',ForgettingFactor=0.99, ...
                YLimits=[-30 10],ShowLegend=true);
    spectrumScope([x_after_adc]);

    % spectrumScope([x_after_subadc(1:nn*2), x_after_adc]);