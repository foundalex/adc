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
       'timeSkew8', '0.8', ...
       'Gain2', '1.2', ...
       'Gain3', '1.3', ...
       'Gain4', '1.4', ...
       'Gain5', '1.5', ...
       'Gain6', '1.6', ...
       'Gain7', '1.7', ...
       'Gain8', '1.8');
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

case 'SNR'
   snr_str = get(curr_obj,'String');
   try
      snr = eval(snr_str);
      state = setfield(state, 'SNR', snr_str);
   catch
      set(curr_obj,'String', state.SNR);
      errordlg('SNR value not a number','Invalid input', 'modal');
   end
   
otherwise
   
end

