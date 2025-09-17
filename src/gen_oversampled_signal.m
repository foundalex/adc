function [s_to_subadc, adc_input, adc_input_int, s_after_subadc, sim_options] = gen_oversampled_signal(sim_options)

    dt = 1/sim_options.Fs;                                                                          % seconds per sample
    t = 0:dt:sim_options.StopTime;                                                                  % seconds
                                                                                                    % step frequency input signal
    sim_options.freq = sim_options.freq + sim_options.step;                                         % frequency of fundamental tone
    sim_options.Z = ceil(sim_options.freq/(sim_options.Fs/sim_options.Inter/2/sim_options.M));      % Nyquist zone

    % Create main signal with noise in double
    s = (cos(2*pi*sim_options.freq*t));
    
    % add noise
    s = awgn(s,sim_options.SNR, "measured");

    % % s + noise integer
    s_fi = fi(s,1,12,11);
    s_int = int16(round(s_fi*2^11));


    % oversampled signal transfer to sub-adc
    offset = 200;
    for i = 1:sim_options.M
        sig = s(i*sim_options.Inter+1+offset:sim_options.M*sim_options.Inter:end);
        sig_int = s_int(i*sim_options.Inter+1+offset:sim_options.M*sim_options.Inter:end);
        if i > 1
            indexx(i-1) = (i*sim_options.Inter*+1+offset);
        end

        adc_input(1:length(sig),i) = sig;
        adc_input_int(1:length(sig_int),i) = sig_int;
    end
    adc_input(length(sig):end,:) = []; 
    
    % adc_input(:,1) = awgn(adc_input(:,1),sim_options.SNR, "measured");
    % adc_input(:,2) = awgn(adc_input(:,2), 0 , "measured");


    adc_input_int(length(sig_int):end,:) = [];

	% исходный сигнал до искажений
	s_to_subadc = zeros(sim_options.M*length(adc_input(:,1)),1);
    s_to_subadc_int = zeros(sim_options.M*length(adc_input_int(:,1)),1);
	for i = 1:sim_options.M
		s_to_subadc(i:sim_options.M:end) = adc_input(:,i); 
        s_to_subadc_int(i:sim_options.M:end) = adc_input_int(:,i);
	end

    %% Add time skew error, gain error
    if sim_options.MODEL_ERROR == true
        % time skew model
        for i = 1:sim_options.M-1
            adc_input_skew = time_skew_func(sim_options.time_skew_array(i), s, indexx(i), sim_options.Inter, sim_options.M); 
            adc_input(1:length(adc_input_skew),i+1) = adc_input_skew;

            adc_input_skew_int = time_skew_func(sim_options.time_skew_array(i), s_int, indexx(i), sim_options.Inter, sim_options.M); 
            adc_input_int(1:length(adc_input_skew_int),i+1) = adc_input_skew_int;
        end
        adc_input(length(adc_input_skew):end,:) = [];
        adc_input_int(length(adc_input_skew_int):end,:) = [];

	    % % model offset error
        % for i = 1:sim_options.M-1
        %     adc_input(:,i+1) = adc_input(:,i+1) + offset_error_array(i);
        % end

        % model gain error
        for i = 1:sim_options.M-1
            adc_input(:,i+1) = adc_input(:,i+1) * sim_options.gain_error_array(i);
            adc_input_int(:,i+1) = int16(fi((adc_input_int(:,i+1)) * sim_options.gain_error_array(i),1,12,0));
        end
    end

    % Main signal with gain error and time skew
	s_after_subadc = zeros(sim_options.M*length(adc_input(:,1)),1);
    s_after_subadc_int = zeros(sim_options.M*length(adc_input_int(:,1)),1);
	for i = 1:sim_options.M
	    s_after_subadc(i:sim_options.M:end) = adc_input(:,i); 
        s_after_subadc_int(i:sim_options.M:end) = adc_input_int(:,i); 
    end

    % save (sprintf(num2str(clock)) + ".mat");
end
%%
% function for model timing skew
function adc_input_skew = time_skew_func(time_skew, s, indexx, Inter, num_adc) 
    adc_input_skew = s(indexx + time_skew*Inter:num_adc*Inter:end);
end