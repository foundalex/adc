function sim_options = ui_read_options

% packet lengths vector, in bits
initial_freq = eval(get(findobj('Tag', 'Initial_Frequency'),'String'))*1000000;

% Number of ADC
num_ADC = eval(get(findobj('Tag', 'Num_ADC'),'String'));  % Num of sub ADC

num_cycles = eval(get(findobj('Tag', 'Number_of_cycles'),'String'));  % Num of cycles modeling

step = eval(get(findobj('Tag', 'step_of_frequency'),'String'))*1000000;  %  Step of frequency

% Channel models
% if get(findobj('Tag', 'AWGN'),'Value')
%    chan_model = 'AWGN';
% elseif get(findobj('Tag', 'ExponentialDecay'),'Value')
%    chan_model = 'ExponentialDecay';
% end
% exp_decay_trms = eval(get(findobj('Tag', 'ExpDecayTrms'),'String'));

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

% Phase Noise
% use_phase_noise = get(findobj('Tag', 'UsePhaseNoise'),'Value');
% phase_noise_dbc = eval(get(findobj('Tag', 'PhaseNoiseDbcLevel'),'String'));
% phase_noise_cfreq = eval(get(findobj('Tag', 'PhaseNoiseCornerFreq'),'String'));
% phase_noise_floor = eval(get(findobj('Tag', 'PhaseNoiseFloor'),'String'));

% Tx Power Spectrum test
% tx_pwr_spectrum_test = get(findobj('Tag', 'TxSpectrumShape'),'Value');

% Synchronization options
% packet_detection = get(findobj('Tag', 'PacketDetection'),'Value');
% fine_time_sync = get(findobj('Tag', 'FineTimeSync'),'Value');
% freq_sync = get(findobj('Tag', 'FreqSync'),'Value');
% pilot_phase_tracking = get(findobj('Tag', 'PilotPhaseTrack'),'Value');
% channel_estimation = get(findobj('Tag', 'ChannelEst'),'Value');
% rx_timing_offset = eval(get(findobj('Tag', 'RxTimingOffset'),'String'));

%Packets per run options
% pkts_per_run = eval(get(findobj('Tag', 'PktsToSimulate'),'String'));


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

