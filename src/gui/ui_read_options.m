function sim_options = ui_read_options

% packet lengths vector, in bits
initial_freq = eval(get(findobj('Tag', 'Initial_Frequency'),'String'))*1000000;

% magnitude = db2mag(eval(get(findobj('Tag', 'magnitude'),'String')));
magnitude = 0;

enable_mask = get(findobj('Tag', 'enable_mask'),'Value');  % enable mask

enable_log = get(findobj('Tag', 'enable_log'),'Value');  % enable log
% Number of ADC
num_ADC = eval(get(findobj('Tag', 'Num_ADC'),'String'));  % Num of sub ADC

num_cycles = eval(get(findobj('Tag', 'Number_of_cycles'),'String'));  % Num of cycles modeling

step = eval(get(findobj('Tag', 'step_of_frequency'),'String'))*1000000;  %  Step of frequency

% Error ADC modeling
error_adc = get(findobj('Tag', 'error_adc'),'Value');


snr = eval(get(findobj('Tag', 'SNR_sub_ADC'),'String'));


% Time Skew Error
for i = 1:8
    time_skew_array(i) = eval(get(findobj('Tag',strcat('timeSkew',string(i))),'String'));
end
% Gain Error
for i = 1:8
    gain_error_array(i) = eval(get(findobj('Tag',strcat('Gain',string(i))),'String'));
end
% Offset Error
for i = 1:8
    offset_error_array(i) = eval(get(findobj('Tag',strcat('Offset',string(i))),'String'));
end


if (enable_log == true)
    diary_name = (['log_file_', num2str(clock), '.txt']);
    diary(diary_name);
end

% Oversampling factor
Inter = 100;
Fs_sub_adc = 1*10^9; % 1 GHz
Fs = Fs_sub_adc * Inter * num_ADC; 


sim_options = struct(                                       ...
   'freq',                          initial_freq,           ...
   'Bit',                           2048,                   ...
   'M',                             num_ADC,                ...
   'num_cycles',                    num_cycles,             ...
   'step',                          step,                   ...
   'SNR',                           snr,                    ...
   'MODEL_ERROR',                   error_adc,              ...
   'time_skew_array',               time_skew_array,        ...
   'gain_error_array',              gain_error_array,       ...
   'offset_error_array',            offset_error_array,     ...
   'Fs',                            Fs,                     ... % Fs all ADC system
   'Fs_sub_adc',                    Fs_sub_adc,             ...
   'StopTime',                      0.001,                  ... % 1 ms need minimum
   'Inter',                         Inter,                  ... % oversampling factor,                      
   'enable_mask',                   enable_mask,            ... % включение наложения маски
   'enable_log',                    enable_log              ... % включение логирования
   );
