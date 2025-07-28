function [yri_cut, te] = fractional_delays(input_signal, M, N_taps, Z)

    n = (1:N_taps);
    ww = blackman(N_taps);
    % % wvtool(blackman(N_taps));
    Nbp = floor(Z/2);
    nn = 1:length(input_signal);

    for i = 1:M-1
        % Method 1
        delay_adc = i/M; % (стр.6,(16)), создаем массив [1:3] из задержанного сигнала ADC0 на различные значения [0.75 0.5 0.25]
        hri = sinc(n - (N_taps - 1) / 2 - delay_adc);
        te(:,i) = hri.' .* ww;
        te(:,i) = te(:,i) ./ sum(te(:,i)); 
        % freqz(te,1);
        yri(:,i) = filter(te(:,i), 1, input_signal); % (стр.6 (15))
        % grpdelay(te(:,i));
        % 

        %Method 2
        hri = designFracDelayFIR(delay_adc, N_taps);
        % Create an FIR filter object
        fdfir = dsp.FIRFilter(hri);
        y = fdfir(input_signal);

        % [gd,f]= grpdelay(hri,1,256,1000000000);  

        % plot_sequences(nn,input_signal, nn,y);
        % plot_sequences(n+delay_adc,input_signal, nn,y);

        % plot([yri(1:100,1), input_signal(1:100)]);

        %% Algorithm for working in different zones of Nyquist
        yhil = hilbert(yri(:,i));

        nn1 = nn + delay_adc;
        cc = cos(2*pi*nn1*Nbp);
        ss = sin(2*pi*nn1*Nbp);

        yri_cos = cc.' .* yri(:,i); 
        yri_sin = ss.' .* imag(yhil); 

        if (mod(Z,2) == 0)
            yric = yri_cos + yri_sin;
        else
            yric = yri_cos - yri_sin;
        end

        yri_cut(:,i) = yric((N_taps-1)/2:end); % удаляем переходные процессы (N-1)/2

    end
    % figure(4);
    % % subplot(2,1,1)
    % % plot([input_signal((N_taps-1)/2:100+35), yri_cut(1:100,1), yri_cut(1:100,2), yri_cut(1:100,3)]);
    % % subplot(2,1,2)
    % plot([input_signal(1:100), yri(1:100,1)]); %, yri(1:100,2), yri(1:100,3)]);
    %     % title('')
    % xlabel('Отсчеты') 
    % ylabel('Амплитуда') 
    % legend('до фильтра дробной задержки','после фильтра дробной задержки')
    % 
    % spectrumScope = spectrumAnalyzer(SampleRate=500000000, ...            
    %             AveragingMethod='exponential',ForgettingFactor=0, ...
    %             YLimits=[-30 10],ShowLegend=true, Method='Welch');
    % 
    % spectrumScope.WindowLength = 2048;
    % spectrumScope.FrequencyResolutionMethod = "window-length";
    % spectrumScope.PlotAsTwoSidedSpectrum=true;
    % spectrumScope.DistortionMeasurements.Enabled = true;
    % 
    % spectrumScope([input_signal(1:4096), yri_cut(1:4096,1)]); %, yri_cut(1:4096,2), yri_cut(1:4096,3)]);
end