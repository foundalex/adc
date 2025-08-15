function sim_options = ui_read_options

% packet lengths vector, in bits
initial_freq = eval(get(findobj('Tag', 'Initial_Frequency'),'String'))*1000000;

% Number of ADC
num_ADC = eval(get(findobj('Tag', 'Num_ADC'),'String'));  % Num of sub ADC

num_cycles = eval(get(findobj('Tag', 'Number_of_cycles'),'String'));  % Num of cycles modeling

step = eval(get(findobj('Tag', 'step_of_frequency'),'String'))*1000000;  %  Step of frequency

% Signal to Noise Rations
snr = eval(get(findobj('Tag', 'SNR'),'String'));

% Error ADC modeling
error_adc = get(findobj('Tag', 'error_adc'),'Value');

% Time Skew Error
for i = 2:8
    time_skew_array(i-1) = eval(get(findobj('Tag',strcat('timeSkew',string(i))),'String'));
end
% Gain Error
for i = 2:8
    gain_error_array(i-1) = eval(get(findobj('Tag',strcat('Gain',string(i))),'String'));
end


Inter = 20;
Fs = 1000000000 * Inter * num_ADC, 

sim_options = struct('freq', initial_freq, ...
   'M', num_ADC, ...
   'num_cycles', num_cycles, ...
   'step', step, ...
   'SNR', snr, ...
   'MODEL_ERROR', error_adc, ...
   'time_skew_array', time_skew_array, ...
   'gain_error_array', gain_error_array, ...
   'Fs',  Fs,           ...    % Fs all ADC system
   'N',  73,            ...    % Number taps filters 
   'StopTime', 0.00001, ...    % seconds
   'Inter', Inter       ...    % oversampling factor
   );

