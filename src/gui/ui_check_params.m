function ui_check_params

persistent state;

if isempty(state)
   state = struct('Num_ADC', '8', ...
       'Initial_Frequency', '20', ...
       'SNR', '65', ...
       'timeSkew2', '0.2', ...
       'timeSkew3', '0.3', ...
       'timeSkew4', '0.4', ...
       'timeSkew5', '0.5', ...
       'timeSkew6', '0.6', ...
       'timeSkew7', '0.7', ...
       'timeSkew8', '0.8');
end

curr_obj = gcbo;
obj_tag = get(curr_obj,'Tag');

% performs logical check on input parameters
switch (obj_tag)
  
   
case 'timeSkew2'
   timeSkew_str = get(curr_obj,'String');
   try
      time_skew2 = eval(timeSkew_str);
      if time_skew2 > 2 | time_skew2 < -2 
         set(curr_obj, 'String', state.timeSkew2);
         errordlg('Time skew 2 не может быть больше 1 или меньше -1 ','Invalid input', 'modal');
      else
         state = setfield(state, 'timeSkew2', timeSkew_str);
      end
   catch
      set(curr_obj,'String', state.timeSkew2);
      errordlg('Time skew value not a valid number','Invalid input', 'modal');
   end
%%
case 'timeSkew3'
   timeSkew_str = get(curr_obj,'String');
   try
      time_skew3 = eval(timeSkew_str);
      if time_skew3 > 2 | time_skew3 < -2 
         set(curr_obj, 'String', state.timeSkew3);
         errordlg('Time skew 3 не может быть больше 1 или меньше -1 ','Invalid input', 'modal');
      else
         state = setfield(state, 'timeSkew3', timeSkew_str);
      end
   catch
      set(curr_obj,'String', state.timeSkew3);
      errordlg('Time skew value not a valid number','Invalid input', 'modal');
   end
%%
case 'timeSkew4'
   timeSkew_str = get(curr_obj,'String');
   try
      time_skew4 = eval(timeSkew_str);
      if time_skew4 > 2 | time_skew4 < -2
         set(curr_obj, 'String', state.timeSkew4);
         errordlg('Time skew 4 не может быть больше 1 или меньше -1 ','Invalid input', 'modal');
      else
         state = setfield(state, 'timeSkew4', timeSkew_str);
      end
   catch
      set(curr_obj,'String', state.timeSkew4);
      errordlg('Time skew value not a valid number','Invalid input', 'modal');
   end
%%
case 'timeSkew5'
   timeSkew_str = get(curr_obj,'String');
   try
      time_skew5 = eval(timeSkew_str);
      if time_skew5 > 2 | time_skew5 < -2 
         set(curr_obj, 'String', state.timeSkew5);
         errordlg('Time skew 5 не может быть больше 1 или меньше -1 ','Invalid input', 'modal');
      else
         state = setfield(state, 'timeSkew5', timeSkew_str);
      end
   catch
      set(curr_obj,'String', state.timeSkew5);
      errordlg('Time skew value not a valid number','Invalid input', 'modal');
   end
%%
case 'timeSkew6'
   timeSkew_str = get(curr_obj,'String');
   try
      time_skew6 = eval(timeSkew_str);
      if time_skew6 > 2 | time_skew6 < -2 
         set(curr_obj, 'String', state.timeSkew6);
         errordlg('Time skew 6 не может быть больше 1 или меньше -1 ','Invalid input', 'modal');
      else
         state = setfield(state, 'timeSkew6', timeSkew_str);
      end
   catch
      set(curr_obj,'String', state.timeSkew6);
      errordlg('Time skew value not a valid number','Invalid input', 'modal');
   end
%%
case 'timeSkew7'
   timeSkew_str = get(curr_obj,'String');
   try
      time_skew7 = eval(timeSkew_str);
      if time_skew7 > 2 | time_skew7 < -2 
         set(curr_obj, 'String', state.timeSkew7);
         errordlg('Time skew 7 не может быть больше 1 или меньше -1 ','Invalid input', 'modal');
      else
         state = setfield(state, 'timeSkew7', timeSkew_str);
      end
   catch
      set(curr_obj,'String', state.timeSkew7);
      errordlg('Time skew value not a valid number','Invalid input', 'modal');
   end
%%
case 'timeSkew8'
   timeSkew_str = get(curr_obj,'String');
   try
      time_skew8 = eval(timeSkew_str);
      if time_skew8 > 2 | time_skew8 < -2 
         set(curr_obj, 'String', state.timeSkew8);
         errordlg('Time skew 8 не может быть больше 1 или меньше -1 ','Invalid input', 'modal');
      else
         state = setfield(state, 'timeSkew8', timeSkew_str);
      end
   catch
      set(curr_obj,'String', state.timeSkew8);
      errordlg('Time skew value not a valid number','Invalid input', 'modal');
   end


% case 'FreqError'
%    freq_err_str = get(curr_obj, 'String');
%    try
%       freq_err = eval(freq_err_str);
%       state = setfield(state, 'FreqError', freq_err_str);
%    catch
%       set(curr_obj,'String', state.FreqError);
%       errordlg('Frequency error value not a number','Invalid input', 'modal');
%    end
% case 'AWGN'
%    awgn_val = get(curr_obj, 'Value');
% 
%    if awgn_val == 1
%       set(curr_obj, 'Enable', 'inactive')
%       exp_decay_chan = findobj('Tag', 'ExponentialDecay');
%       set(exp_decay_chan, 'Enable', 'on')
%       set(exp_decay_chan, 'Value', 0);
%    end
% case 'ExponentialDecay'
%    exp_decay_val = get(curr_obj, 'Value');
%    if exp_decay_val == 1      
%       set(curr_obj, 'Enable', 'inactive');
%       awgn_chan = findobj('Tag', 'AWGN');
%       set(awgn_chan, 'Enable', 'on')
%       set(awgn_chan, 'Value', 0);
%    end
% case 'ExpDecayTrms'
%    exp_decay_trms_str = get(curr_obj,'String');
%    try
%       exp_decay_trms = eval(exp_decay_trms_str);
%       if exp_decay_trms < 0
%          set(curr_obj, 'String', state.ExpDecayTrms);
%          errordlg('Exponential decay T rms cannot be negative','Invalid input', 'modal');
%       else
%          state = setfield(state, 'ExpDecayTrms', exp_decay_trms_str);      
%       end
%    catch
%       set(curr_obj,'String', state.ExpDecayTrms);
%       errordlg('Exponential decay T rms value not a number','Invalid input', 'modal');
%    end
case 'SNR'
   snr_str = get(curr_obj,'String');
   try
      snr = eval(snr_str);
      state = setfield(state, 'SNR', snr_str);
   catch
      set(curr_obj,'String', state.SNR);
      errordlg('SNR value not a number','Invalid input', 'modal');
   end
% case 'PhaseNoiseDbcLevel'
%    phase_noise_dbc_str = get(curr_obj,'String');
%    try
%       phase_noise_dbc = eval(phase_noise_dbc_str);
%       if phase_noise_dbc > 0
%          set(curr_obj, 'String', state.PhaseNoisedBc);
%          errordlg('Phase noise dBc level must be negative', 'Invalid input', 'modal');
%       else
%          state = setfield(state, 'PhaseNoisedBc', phase_noise_dbc_str);      
%       end
%    catch
%       set(curr_obj,'String', state.PhaseNoisedBc);
%       errordlg('Phase noise dBc value not a number','Invalid input', 'modal');
%    end
% case 'PhaseNoiseCornerFreq'
%    phase_noise_cfreq_str = get(curr_obj,'String');
%    try
%       phase_noise_cfreq = eval(phase_noise_cfreq_str);
%       if phase_noise_cfreq < 0
%          set(curr_obj, 'String', state.PhaseNoiseCFreq);
%          errordlg('Phase noise corner frequency must be positive','Invalid input', 'modal');
%       else
%          state = setfield(state, 'PhaseNoiseCFreq', phase_noise_cfreq_str);      
%       end
%    catch
%       set(curr_obj,'String', state.PhaseNoiseCFreq);
%       errordlg('Phase noise corner frequency value not a number','Invalid input', 'modal');
%    end
% case 'PhaseNoiseFloor'
%    phase_noise_floor_str = get(curr_obj,'String');
%    try
%       phase_noise_floor = eval(phase_noise_floor_str);
%       if phase_noise_floor > 0
%          set(curr_obj, 'String', state.PhaseNoiseFloor);
%          errordlg('Phase noise floor must be negative','Invalid input', 'modal');
%       else
%          state = setfield(state, 'PhaseNoiseFloor', phase_noise_floor_str);      
%       end
%    catch
%       set(curr_obj,'String', state.PhaseNoiseFloor);
%       errordlg('Phase noise floor level value not a number','Invalid input', 'modal');
%    end
% case 'PacketDetection'
%    pkt_det = get(curr_obj,'Value');
%    if pkt_det == 1
%       fine_time_sync = findobj('Tag', 'FineTimeSync');
%       freq_sync = findobj('Tag', 'FreqSync');
%       pilot_phase_track = findobj('Tag', 'PilotPhaseTrack');
%       channel_est = findobj('Tag', 'ChannelEst');
% 
%       set(fine_time_sync, 'Value', 1);
%       set(freq_sync, 'Value', 1);
%       set(pilot_phase_track, 'Value', 1);
%       set(channel_est, 'Value', 1);
%    end
% case 'FineTimeSync'
%    fine_time_sync = get(curr_obj,'Value');
%    if fine_time_sync == 1
%       freq_sync = findobj('Tag', 'FreqSync');
%       pilot_phase_track = findobj('Tag', 'PilotPhaseTrack');
%       channel_est = findobj('Tag', 'ChannelEst');
% 
%       set(freq_sync, 'Value', 1);
%       set(pilot_phase_track, 'Value', 1);
%       set(channel_est, 'Value', 1);
%    else
%       packet_detection = findobj('Tag', 'PacketDetection');
%       set(packet_detection, 'Value', 0);
%    end
% case 'FreqSync'
%    freq_sync = get(curr_obj,'Value');
%    if freq_sync == 1
%       pilot_phase_track = findobj('Tag', 'PilotPhaseTrack');
%       channel_est = findobj('Tag', 'ChannelEst');
% 
%       set(pilot_phase_track, 'Value', 1);
%       set(channel_est, 'Value', 1);
%    else
%       packet_detection = findobj('Tag', 'PacketDetection');
%       fine_time_sync = findobj('Tag', 'FineTimeSync');
% 
%       set(packet_detection, 'Value', 0);
%       set(fine_time_sync, 'Value', 0);
%    end
% case 'PilotPhaseTrack'
%    pilot_phase_track = get(curr_obj,'Value');
%    if pilot_phase_track == 1
%       channel_est = findobj('Tag', 'ChannelEst');
% 
%       set(channel_est, 'Value', 1);
%    else
%       packet_detection = findobj('Tag', 'PacketDetection');
%       fine_time_sync = findobj('Tag', 'FineTimeSync');
%       freq_sync = findobj('Tag', 'FreqSync');
% 
%       set(packet_detection, 'Value', 0);
%       set(fine_time_sync, 'Value', 0);
%       set(freq_sync, 'Value', 0);
%    end
% case 'ChannelEst'
%    channel_est = get(curr_obj,'Value');
%    if channel_est == 0
%       packet_detection = findobj('Tag', 'PacketDetection');
%       fine_time_sync = findobj('Tag', 'FineTimeSync');
%       freq_sync = findobj('Tag', 'FreqSync');
%       pilot_phase_track = findobj('Tag', 'PilotPhaseTrack');
% 
%       set(packet_detection, 'Value', 0);
%       set(fine_time_sync, 'Value', 0);
%       set(freq_sync, 'Value', 0);
%       set(pilot_phase_track, 'Value', 0);
%    end
% case 'RxTimingOffset'
%    rx_timing_offset_str = get(curr_obj,'String');
%    try
%       rx_timing_offset = eval(rx_timing_offset_str);
%       if rx_timing_offset > 0
%          set(curr_obj,'String', state.RxTimingOffset);
%          errordlg('Rx timing offset positive','Invalid input', 'modal');
%       else
%          state = setfield(state, 'RxTimingOffset', rx_timing_offset_str);
%       end
%    catch
%       set(curr_obj,'String', state.RxTimingOffset);
%       errordlg('Rx timing offset value not a number','Invalid input', 'modal');
%    end
% case 'PktsToSimulate'
%    pkts_to_simulate_str = get(curr_obj, 'String');
%    try
%       pkts_to_simulate = eval(pkts_to_simulate_str);
%       if pkts_to_simulate < 0
%          set(curr_obj,'String', state.PktsPerRun);
%          errordlg('Packets to simulate cannot be negative','Invalid input', 'modal');
%       else
%          state = setfield(state, 'PktsPerRun', pkts_to_simulate_str);
%       end
%    catch
%       set(curr_obj,'String', state.PktsPerRun);
%       errordlg('Packets to simulate value not a number','Invalid input', 'modal');
%    end
otherwise
   
end

