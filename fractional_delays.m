function [yri_cut, te] = fractional_delays(input_signal, M, N_taps, Z)

    n = (1:N_taps);

    ww = blackman(N_taps);
    % % wvtool(blackman(N_taps));
    Nbp = floor(Z/2);
    nn = 1:length(input_signal);

    del_proc = ((N_taps-1)/2)+1;

    for i = 1:M-1
        % Method 1
        delay_adc = i/M; % (стр.6,(16)), создаем массив [1:3] из задержанного сигнала ADC0 на различные значения [0.75 0.5 0.25]

        D = ((N_taps-1)/2) + delay_adc;
        hri = sinc(n - ((N_taps - 1) / 2) - delay_adc);
        hri_m = sin(pi * (n - D)) ./ (pi * (n - D));
        mean(grpdelay(hri_m,1,256,1))
        grpdelay(hri,1,256,1);

        hri_m = hri_m.' .* ww;
        hri_m = hri_m ./ sum(hri_m); 

        te(:,i) = hri.' .* ww;
        te(:,i) = te(:,i) ./ sum(te(:,i)); 

        % plot([te(:,i), hri_m]);
        % freqz(te,1);
        % freqz(hri_m,1);
        
        yri(:,i) = filter(te(:,i), 1, input_signal); % (стр.6 (15))
        yri1(:,i) = filter(hri_m, 1, input_signal); % (стр.6 (15))
    
        plot([[zeros(35,1); input_signal(1:30)], yri(1:30+35,1), yri1(1:30+35,1)]);

        %Method 2
        % hri = designFracDelayFIR(delay_adc, N_taps);
        % % Create an FIR filter object
        % fdfir = dsp.FIRFilter(hri);
        % y(:,i) = fdfir(input_signal);

     

        %% Algorithm for working in different zones of Nyquist
        yhil = hilbert(yri(:,i));
        yhil_imag(:,i) = imag(yhil);
        nn1 = nn + delay_adc;

        cc = cos(2*pi*nn1*Nbp);
        ss = sin(2*pi*nn1*Nbp);


        yri_cos = cc.' .* yri(:,i); 
        yri_sin = ss.' .* yhil_imag(:,i); 



        if (mod(Z,2) == 0)
            yric = yri_cos + yri_sin;
        else
            yric = yri_cos - yri_sin;
        end

        yri_cut(:,i+1) = yri1(1:end,i); % удаляем переходные процессы (N-1)/2
    end


    a1 = imag(hilbert(input_signal));
    yri_cut(:,1) = input_signal;

    plot([[zeros(37,1); input_signal(1:23)], [0; yri1(1:59,1)], yri1(1:60,2)]); %, yri1(1:60,3)]);
    xlabel('Отсчеты') 
    ylabel('Амплитуда') 
    legend('до фильтра дробной задержки','задержка на 0.25', 'задержка на 0.5', 'задержка на 0.75');

    % figure(4);
    % % subplot(2,1,1)
    % plot([input_signal(1:100), yri_cut(1:100,1), yri_cut(1:100,2), yri_cut(1:100,3)]);
    % subplot(2,1,2)
    % plot([input_signal(1:100), yri(1:100,1), yri(1:100,2), yri(1:100,3)]);
    %     % title('')
    % xlabel('Отсчеты') 
    % ylabel('Амплитуда') 
    % legend('до фильтра дробной задержки','после фильтра дробной задержки')
    % 
    % spectrumScope = spectrumAnalyzer(SampleRate=1000000000, ...            
    %             AveragingMethod='exponential',ForgettingFactor=0, ...
    %             YLimits=[-30 10],ShowLegend=true, Method='Welch');
    % % 
    % spectrumScope.WindowLength = 2048;
    % spectrumScope.FrequencyResolutionMethod = "window-length";
    % spectrumScope.PlotAsTwoSidedSpectrum=true;
    % spectrumScope.DistortionMeasurements.Enabled = true;
    % 
    % spectrumScope([input_signal(1:4096), yri_cut(1:4096,1), yri_cut(1:4096,2), yri_cut(1:4096,3)]);
end