function [s_to_subadc, s_to_subadc_int, adc_input, adc_input_int, s_after_subadc, s_after_subadc_int] = gen_oversampled_signal(M, Fs, freq, SNR, Inter, StopTime, ...
    MODEL_ERROR, time_skew_array, gain_error_array)

    dt = 1/Fs;                                                                          % seconds per sample
    t = 0:dt:StopTime;                                                                  % seconds
                                                                                        % step frequency input signal
    % Create main signal with noise in double
    s = (cos(2*pi*freq*t));
    
    % add noise
    s = awgn(s, SNR, "measured");

    % % s + noise integer
    s_fi = fi(s,1,12,11);
    s_int = int16(round(s_fi*2^11));

    % s_int = int16(floor(s*2^11));


    % oversampled signal transfer to sub-adc
    offset = 0;
    for i = 1:M
        sig = s(i*Inter+1+offset:M*Inter:end); % 10, 20, 30
        sig_int = s_int(i*Inter+1+offset:M*Inter:end);
        if i > 1
            indexx(i-1) = (i*Inter*+1+offset);
        end

        adc_input(1:length(sig),i) = sig;
        adc_input_int(1:length(sig_int),i) = sig_int;
    end
    adc_input(length(sig):end,:) = []; 
    
    % adc_input(:,1) = awgn(adc_input(:,1),sim_options.SNR, "measured");
    % adc_input(:,2) = awgn(adc_input(:,2), 0 , "measured");


    adc_input_int(length(sig_int):end,:) = [];

    % figure(2);
    % subplot(2,1,1)
    % plot(s(1:200));
    % subplot(2,1,2)
    % plot(adc_input_int(1:200,1));
    % % 
    % spectrumScope = spectrumAnalyzer(SampleRate=Fs/Inter/4, ...            
    %         AveragingMethod='exponential',ForgettingFactor=0.99, ...
    %         YLimits=[-30 10],ShowLegend=true);
    % 
    % spectrumScope([adc_input(:,1)]);


	% исходный сигнал до искажений
	s_to_subadc = zeros(M*length(adc_input(:,1)),1);
    s_to_subadc_int = zeros(M*length(adc_input_int(:,1)),1);
	for i = 1:M
		s_to_subadc(i:M:end) = adc_input(:,i); 
        s_to_subadc_int(i:M:end) = adc_input_int(:,i);
	end

    %% Add time skew error, gain error
    if MODEL_ERROR == true
        % time skew model
        for i = 1:M-1
            adc_input_skew = time_skew_func(time_skew_array(i), s, indexx(i), Inter, M); 
            adc_input(1:length(adc_input_skew),i+1) = adc_input_skew;

            adc_input_skew_int = time_skew_func(time_skew_array(i), s_int, indexx(i), Inter, M); 
            adc_input_int(1:length(adc_input_skew_int),i+1) = adc_input_skew_int;
        end
        adc_input(length(adc_input_skew):end,:) = [];
        adc_input_int(length(adc_input_skew_int):end,:) = [];

	    % % model offset error
        % for i = 1:sim_options.M-1
        %     adc_input(:,i+1) = adc_input(:,i+1) + offset_error_array(i);
        % end

        % model gain error
        for i = 1:M-1
            adc_input(:,i+1) = adc_input(:,i+1) * gain_error_array(i);
            adc_input_int(:,i+1) = int16(fi((adc_input_int(:,i+1)) * gain_error_array(i),1,12,0));
        end
    end

    % Main signal with gain error and time skew
	s_after_subadc = zeros(M*length(adc_input(:,1)),1);
    s_after_subadc_int = zeros(M*length(adc_input_int(:,1)),1);
	for i = 1:M
	    s_after_subadc(i:M:end) = adc_input(:,i); 
        s_after_subadc_int(i:M:end) = adc_input_int(:,i); 
    end


    % save (sprintf(num2str(clock) + ".mat"));
    % load ('2025             10              8             13             24         10.929.mat'); % 70 SNR 70 MHz
    % load ('2025             10             17             15             47          4.423.mat'); % 60 SNR 70 MHz

    % load ('2025             10             29             13             27          1.297.mat'); % 0 SNR 298 MHz
    % load ('2025             10             28             13             31         30.067.mat'); % 0 SNR 298 MHz

    % load ('2025             10             30             12
    % 35          34.93.mat'); % 70 SNR 741 MHz 999 samples

     load ('2025             11              1             15             44         59.704.mat'); % 70 SNR 741 MHz 9999 samples
end

%%
% function for model timing skew
function adc_input_skew = time_skew_func(time_skew, s, indexx, Inter, num_adc) 
    adc_input_skew = s(indexx + time_skew*Inter:num_adc*Inter:end);
end