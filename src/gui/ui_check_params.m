function ui_check_params

persistent state;

if isempty(state)
   state = struct('Num_ADC', '4', ...
       'Initial_Frequency', '50', ...
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
       'Gain5', '1.2', ...
       'Gain6', '1.3', ...
       'Gain7', '1.4', ...
       'Gain8', '0.9');
end

curr_obj = gcbo;
obj_tag = get(curr_obj,'Tag');

% performs logical check on input parameters
switch (obj_tag)

  %% проверка входной частоты сигнала. Алгоритм калибровки работает с сигналом >= 50 МГц
  case 'Initial_Frequency'
  init_freq_str = get(curr_obj,'String');
   try
      init_freq = eval(init_freq_str);
      if init_freq < 50
         set(curr_obj, 'String', state.Initial_Frequency);
         errordlg('Частота сигнала не может быть меньше 50 МГц','Invalid input', 'modal');
      else
         state = setfield(state, 'Initial_Frequency', init_freq_str);
      end
   catch
      set(curr_obj,'String', state.Initial_Frequency);
      errordlg('Initial frequency value not a valid number','Invalid input', 'modal');
   end
 

  %% проверка ошибки gain. Алгоритм калибровки работает с 0.1 =< gain error =< 1.4 

   case 'Gain2'
   Gain2_str = get(curr_obj,'String');
   try
      Gain2i = eval(Gain2_str);
      if Gain2i > 1.4 | Gain2i < 0.1 
         set(curr_obj, 'String', state.Gain2);
         errordlg('Gain ADC 2 не может быть больше 1.4 или меньше 0.1', 'Invalid input', 'modal');
      else
         state = setfield(state, 'Gain2', Gain2_str);
      end
   catch
      set(curr_obj,'String', state.Gain2);
      errordlg('Gain 2 ADC value not a valid number','Invalid input', 'modal');
   end


   case 'Gain3'
   Gain3_str = get(curr_obj,'String');
   try
      Gain3i = eval(Gain3_str);
      if Gain3i > 1.4 | Gain3i < 0.1 
         set(curr_obj, 'String', state.Gain3);
         errordlg('Gain ADC 3 не может быть больше 1.4 или меньше 0.1', 'Invalid input', 'modal');
      else
         state = setfield(state, 'Gain3', Gain3_str);
      end
   catch
      set(curr_obj,'String', state.Gain3);
      errordlg('Gain 3 ADC value not a valid number','Invalid input', 'modal');
   end


   case 'Gain4'
   Gain4_str = get(curr_obj,'String');
   try
      Gain4i = eval(Gain4_str);
      if Gain4i > 1.4 | Gain4i < 0.1 
         set(curr_obj, 'String', state.Gain4);
         errordlg('Gain ADC 4 не может быть больше 1.4 или меньше 0.1', 'Invalid input', 'modal');
      else
         state = setfield(state, 'Gain4', Gain4_str);
      end
   catch
      set(curr_obj,'String', state.Gain4);
      errordlg('Gain 4 ADC value not a valid number','Invalid input', 'modal');
   end


   case 'Gain5'
   Gain5_str = get(curr_obj,'String');
   try
      Gain5i = eval(Gain5_str);
      if Gain5i > 1.4 | Gain5i < 0.1 
         set(curr_obj, 'String', state.Gain5);
         errordlg('Gain ADC 5 не может быть больше 1.4 или меньше 0.1', 'Invalid input', 'modal');
      else
         state = setfield(state, 'Gain5', Gain5_str);
      end
   catch
      set(curr_obj,'String', state.Gain5);
      errordlg('Gain 5 ADC value not a valid number','Invalid input', 'modal');
   end


   case 'Gain6'
   Gain6_str = get(curr_obj,'String');
   try
      Gain6i = eval(Gain6_str);
      if Gain6i > 1.4 | Gain6i < 0.1 
         set(curr_obj, 'String', state.Gain6);
         errordlg('Gain ADC 6 не может быть больше 1.4 или меньше 0.1', 'Invalid input', 'modal');
      else
         state = setfield(state, 'Gain6', Gain6_str);
      end
   catch
      set(curr_obj,'String', state.Gain6);
      errordlg('Gain 6 ADC value not a valid number','Invalid input', 'modal');
   end


   case 'Gain7'
   Gain7_str = get(curr_obj,'String');
   try
      Gain7i = eval(Gain7_str);
      if Gain7i > 1.4 | Gain7i < 0.1 
         set(curr_obj, 'String', state.Gain7);
         errordlg('Gain ADC 7 не может быть больше 1.4 или меньше 0.1', 'Invalid input', 'modal');
      else
         state = setfield(state, 'Gain7', Gain7_str);
      end
   catch
      set(curr_obj,'String', state.Gain7);
      errordlg('Gain 7 ADC value not a valid number','Invalid input', 'modal');
   end


   case 'Gain8'
   Gain8_str = get(curr_obj,'String');
   try
      Gain8i = eval(Gain8_str);
      if Gain8i > 1.4 | Gain8i < 0.1 
         set(curr_obj, 'String', state.Gain8);
         errordlg('Gain ADC 8 не может быть больше 1.4 или меньше 0.1', 'Invalid input', 'modal');
      else
         state = setfield(state, 'Gain8', Gain8_str);
      end
   catch
      set(curr_obj,'String', state.Gain8);
      errordlg('Gain 8 ADC value not a valid number','Invalid input', 'modal');
   end



%%

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

