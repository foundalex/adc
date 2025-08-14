function [x_to_subadc, adc_input, x_after_subadc, sim_options] = gen_oversampled_signal(sim_options)

    dt = 1/sim_options.Fs;                                                  % seconds per sample
    t = 0:dt:sim_options.StopTime;                                          % seconds
                                                                            % step frequency input signal
    sim_options.freq = sim_options.freq + sim_options.step;                             % frequency of fundamental tone
    sim_options.Z = ceil(sim_options.freq/(sim_options.Fs/sim_options.Inter/2/sim_options.M));      % Nyquist zone

    % Create main signal with noise
    s = 0.75*cos(2*pi*sim_options.freq*t);
    noise = awgn(s,sim_options.SNR);
    s = s + noise;

    % oversampled signal transfer to sub-adc
    offset = 200;
    for i = 1:sim_options.M
        sig = s(i*sim_options.Inter+1+offset:sim_options.M*sim_options.Inter:end);
        if i > 1
            indexx(i-1) = (i*sim_options.Inter*+1+offset);
        end

        adc_input(1:length(sig),i) = sig;
    end
    adc_input(length(sig):end,:) = []; 

	% исходный сигнал до искажений
	x_to_subadc = zeros(sim_options.M*length(adc_input(:,1)),1);
	for i = 1:sim_options.M
		x_to_subadc(i:sim_options.M:end) = adc_input(:,i); 
	end

    %% Add time skew error, gain error
    if sim_options.MODEL_ERROR == true
        % time skew model
        for i = 1:sim_options.M-1
            adc_input_skew = time_skew_func(sim_options.time_skew_array(i), s, indexx(i), sim_options.Inter, sim_options.M); 
            adc_input(1:length(adc_input_skew),i+1) = adc_input_skew;
        end
        adc_input(length(adc_input_skew):end,:) = [];

        % figure(2);
        % plot([adc_input(1:50,1), adc_input(1:50,2)]);

	    % % model offset error
        % for i = 1:M-1
        %     adc_input(:,i+1) = adc_input(:,i+1) + offset_error_array(i);
        % end

        % model gain error
        for i = 1:sim_options.M-1
            adc_input(:,i+1) = adc_input(:,i+1) * sim_options.gain_error_array(i);
        end
    end

    % Main signal with gain error and time skew
	x_after_subadc = zeros(sim_options.M*length(adc_input(:,1)),1);
	for i = 1:sim_options.M
	    x_after_subadc(i:sim_options.M:end) = adc_input(:,i); 
    end

end

%%
% function for model timing skew
function adc_input_skew = time_skew_func(time_skew, s, indexx, Inter, num_adc) 
    adc_input_skew = s(indexx + time_skew*Inter:num_adc*Inter:end);
end