clear all;

FrameSize = 10;
NIter = 100;

rls1 = dsp.RLSFilter('Length',10);
    % 'Method','Householder RLS');
filt = dsp.FIRFilter('Numerator',...
    fir1(10,[.5,.75]));

sinewave = dsp.SineWave('Frequency',0.01,...
    'SampleRate',1,...
    'SamplesPerFrame',FrameSize);
% scope = timescope('LayoutDimensions',[2 1],...
%     'NumInputPorts',2, ...
%     'TimeUnits','Seconds',...
%     'YLimits',[-2.5 2.5], ...
%     'BufferLength',2*FrameSize*NIter,...
%     'ActiveDisplay',1,...
%     'ShowLegend',true,...
%     'ChannelNames',{'Noisy signal'},...
%     'ActiveDisplay',2,...
%     'ShowLegend',true,...
%     'ChannelNames',{'Error signal'});

xx = zeros(1000,1);
dd = zeros(1000,1);
yy = zeros(1000,1);

%%
% parameter estimator menggunakan least square

% Y = ([1 : 100] + rand(1, 1))'; % proses markov orde 1
% X = ([100 : 199] + rand(1, 1))';
% st_time = tic();
% % jumlah parameter
% np = 2;
% % orde model
% nmdl = np / 2;
% phi = [];
% N = 40;
% for m = nmdl + 1 : N - nmdl
%     phi(m, :) = [-Y(m-1:-1:m-nmdl, :)' X(m-1:-1:m-nmdl, :)'];
% end
% Y2 = Y(nmdl + 1 : N);
% % Parameter yang diestimasi
% beta = inv(phi'*phi) * phi' * Y2;
% Y3 = phi * beta;
% err = Y2 - Y3;
% emse = mse(err);
% t = 1 : length(Y3);
% figure1;
% plot(t, Y2, t, Y3);
% title('Data Asli vs Estimasi');
% figure2;
% plot(t, err);
% title('Galat');
% cnsmd_time = toc(st_time);



for k = 1:NIter
    x = randn(FrameSize,1);
    d = filt(x) + sinewave();
    [y,e] = rls1(x,d);
    w = rls1.Coefficients;

    % scope(d,e)
    %%
    % Запись полного сигнала
    xx((k-1)*length(x)+1:k*length(x)) = x;
    dd((k-1)*length(d)+1:k*length(d)) = d;
    yy((k-1)*length(y)+1:k*length(y)) = y;
end


function [n,x,s,fs] = genData(numPts, freq, filt, nVar, SNR)
    % Generate time values
    t = linspace(0,1,numPts)';
    fs = numPts;
    
    % Generate tone
    s = sin(2*pi*freq*t);
    
    % Generate noise
    n = sqrt(nVar)*randn(numPts,1);
    
    % Filter noise
    addnoise = filter(filt, 1, n);
    
    % Plot filter
    freqz(filt,1,1000)
    
    % Adjust SNR of tone
    s = s/sqrt(var(s)/(10^(SNR/10)*var(n)));
    disp(['Calculated SNR = ' num2str(10*log10(var(s)/var(n)))])
    
    % Add noise to signal
    x = s + addnoise;
    
end

numPts  = 1000;             % number of points to generate
freq    = 100;              % frequency of fundamental tone
filtord = 4;                % filter order
filt    = rand(filtord, 1); % filter coefficients
nVar    = 1;                % white noise variance
SNR     = -20;              % signal to noise ratio of tone
    
% Generate the data!
[n,x,s,fs] = genData(numPts, freq, filt, nVar, SNR);

rr = (1:50);
yy1 = (rr*10).';
x3 = zeros(10,10);


x3(1,:) = x(1:10);
x3(:,1) = x(1:10).';

w1 = inv(x3'*x3) * x3' * s(1:10);
st = lsqr(x3, s(1:10));


y_out = filter(w1.', 1, x(1:10));


for j = 11:length(s)
    x3(1,:) = circshift(x3(1,:),1);
    
    for i = 2:10
        x3(i,:) = circshift(x3(i-1,:),-1);
    end
    x3(1,1) = x(j);

    w1 = inv(x3'*x3) * x3' * s(j-10:j-1);
    % st = lsqr(x3, s(j-10:j-1));

    y_out = filter(w1.', 1, x(j-10:j-1));
    y_out1((j-10:j-1)) = y_out;
   
end

 y_out1 =  y_out1.'

plot([s(1:999),x(1:999), y_out1])
% title('Отношение между отсчетами I-составляющей')
xlabel('Номер отсчета') 
ylabel('Амплитуда') 
legend({'Исходный сигнал','Сигнал + Шум','Сигнал с выхода адаптивного фильтра'},'Location','northeast')
%%
% subplot(2,1,1)
% plot([xx, dd, yy]);
% subplot(2,1,2)
% plot(e)



% release(scope)

    %--------------------------------------------------------------------------
% Filtering
%--------------------------------------------------------------------------
% Filter Parameters
p       = 1;                % filter order
lambda  = 1.0;              % forgetting factor
laminv  = 1/lambda;
delta   = 1.0;              % initialization parameter
% Filter Initialization
w       = zeros(p,1);       % filter coefficients
P       = delta*eye(p);     % inverse correlation matrix
e       = x*0;              % error signal


x2 = zeros(4,4);
y1 = zeros(1,4);

for m = p:length(x)
    % Acquire chunk of data
    y = n(m:-1:m-p+1);
    % Error signal equation
    e(m) = x(m)-w'*y;
    r = x(m);
    % Parameters for efficiency
    Pi = P*y;
    
    % Filter gain vector update
    k = (Pi)/(lambda+y'*Pi);
    % Inverse correlation matrix update
    P = (P - k*y'*P)*laminv;
    % Filter coefficients adaption
    w = w + k*e(m);
    % Counter to show filter is working
    %if mod(m,1000) == 0
    %    disp([num2str(m/1000) ' of ' num2str(floor(length(x)/1000))])
    %end

    % x1 = x(m:-1:m-p+1);
    x1 = x(1:4);
    y1 = n(1:4);
    x3 = zeros(2,2);
    y3 = zeros(2,1);
    x33 = zeros(2,2);
    y33 = zeros(2,1);
    % for i = 1:4
    %     % x1 = [x(i), x1(1,1), x1(1,2), x1(1,3); x1(1,1), x1(2,1), x1(2,2), x1(2,3); ...
    %     %     x1(2,1), x1(3,1), x1(3,2), x1(3,3);...
    %     %     x1(3,1), x1(4,1), x1(4,2), x1(4,3)];
    %     if i == 1
    %         x2 = [x1(1); 0; 0; 0];
    %     elseif i == 2
    %         x2 = [x1(1); x1(2); 0; 0];
    %     elseif i == 3
    %         x2 = [x1(1); x1(2); x1(3); 0];
    %     elseif i == 4
    %         x2 = [x1(1); x1(2); x1(3); x1(4)];
    %     end
    % 
    %     if i == 1
    %         y1 = [n(1); 0; 0; 0];
    %     elseif i == 2
    %         y1 = [n(1); n(2); 0; 0];
    %     elseif i == 3
    %         y1 = [n(1); n(2); n(3); 0];
    %     elseif i == 4
    %         y1 = [n(1); n(2); n(3); n(4)];
    %     end

       % y1 = [n(i); y1(1); y1(2); y1(3)];
     
       % a = lsqr(x1(1), y1(1));
       % b = lsqr(x1(2), y1(2));

       % shift from first order
       % for i = 1:2
       %  x33 = [x(i), x33(1); x33(1,1), x33(2,1)];
       %  y33 = [n(i), y33(1); y33(1,1), y33(2,1)];
       %  w2 = inv(x33'*x33) * x33' * y33;
       % end


       % work
       for i = 1:2
        x3 = [x(i), x(i+1); x(i+1), x(i+2)];
        % x3 = [x(i), x(2); x(2), 0];
        % y3 = [n(1); n(2)];
        y3 = [n(i); n(i+1)];
        w1 = inv(x3'*x3) * x3' * y3;
        st = lsqr(x3, y3);
       end
 
 end
