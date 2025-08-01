function yri_cut = fractional_delays(input_signal, M, N_taps, Z)

    ww = blackman(N_taps);
    Nbp = floor(Z/2);
    nn = 1:length(input_signal);
    n = (0:1:N_taps-1);
    del_proc = ((N_taps-1)/2);

    for i = 1:M-1
        % Method 1
        delay_adc = i/M; % (стр.6,(16)), создаем массив [1:3] из задержанного сигнала ADC0 на различные значения [0.75 0.5 0.25]
        D = ((N_taps-1)/2) - delay_adc;
        % hri = sinc(n - ((N_taps - 1) / 2) - delay_adc);
        hri_m(:,i) = sin(pi * (n - D)) ./ (pi * (n - D));
        % grpdelay(hri_m(:,i),1,256,1);
        % mean(grpdelay(hri_m(:,i),1,256,1))

        hri_m(:,i) = hri_m(:,i) .* ww;
        hri_m(:,i) = hri_m(:,i) ./ sum(hri_m(:,i)); 

        % freqz(hri_m(:,i),1);   
        yri(:,i) = filter(hri_m(:,i), 1, input_signal); % (стр.6 (15))

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
        yri_cut(:,i+1) = yric(((N_taps-1)/2)+1:end);

    
    end
    yri_cut(:,1) = input_signal(1:end - ((N_taps-1)/2));

    % plot([[zeros(36,1); input_signal(1:34)], yri(1:70,1), yri(1:70,2), yri(1:70,3)]);

    % spectrumScope = spectrumAnalyzer(SampleRate=1000000000, ...            
    %             AveragingMethod='exponential',ForgettingFactor=0, ...
    %             YLimits=[-30 10],ShowLegend=true, Method='Welch');
    % spectrumScope.WindowLength = 2048;
    % spectrumScope.FrequencyResolutionMethod = "window-length";
    % spectrumScope.PlotAsTwoSidedSpectrum=true;
    % spectrumScope.DistortionMeasurements.Enabled = true;
    % 
    % spectrumScope([input_signal(100:4096+100,1), input_signal(100:4096+100,2), yri_cut(100:4096+100,2)]);
end