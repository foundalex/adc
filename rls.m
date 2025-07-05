FrameSize = 10;
NIter = 100;

rls1 = dsp.RLSFilter('Length',11);
    % 'Method','Householder RLS');
filt = dsp.FIRFilter('Numerator',...
    fir1(10,[.5,.75]));

sinewave = dsp.SineWave('Frequency',0.01,...
    'SampleRate',1,...
    'SamplesPerFrame',FrameSize);
scope = timescope('LayoutDimensions',[2 1],...
    'NumInputPorts',2, ...
    'TimeUnits','Seconds',...
    'YLimits',[-2.5 2.5], ...
    'BufferLength',2*FrameSize*NIter,...
    'ActiveDisplay',1,...
    'ShowLegend',true,...
    'ChannelNames',{'Noisy signal'},...
    'ActiveDisplay',2,...
    'ShowLegend',true,...
    'ChannelNames',{'Error signal'});

xx = zeros(1000,1);
dd = zeros(1000,1);
yy = zeros(1000,1);

%% пример из книги
x23 = [2 1 0; 1 2 0.1];
d3 = [1; -1; 0];

% x23*x23(transponse)
x23_m = x23*x23.';
% x23 * d3
x23_m_d3 = x23*d3;
ww = inv(x23_m) * x23_m_d3;

%%
% parameter estimator menggunakan least square

Y = ([1 : 100] + rand(1, 1))'; % proses markov orde 1
X = ([100 : 199] + rand(1, 1))';
st_time = tic();
% jumlah parameter
np = 2;
% orde model
nmdl = np / 2;
phi = [];
N = 40;
for m = nmdl + 1 : N - nmdl
    phi(m, :) = [-Y(m-1:-1:m-nmdl, :)' X(m-1:-1:m-nmdl, :)'];
end
Y2 = Y(nmdl + 1 : N);
% Parameter yang diestimasi
beta = inv(phi'*phi) * phi' * Y2;
Y3 = phi * beta;
err = Y2 - Y3;
emse = mse(err);
t = 1 : length(Y3);
figure;
plot(t, Y2, t, Y3);
title('Data Asli vs Estimasi');
figure;
plot(t, err);
title('Galat');
cnsmd_time = toc(st_time);


%%
x1 = [];

for k = 1:NIter
    x = randn(FrameSize,1);
    d = filt(x) + sinewave();
    [y,e] = rls1(x,d);
    w = rls1.Coefficients;
    scope(d,e)

    %%
    x1(:,k) = x.'; 
    w1 = inv(x1'*x1) * x1' * d;
    err1 = d - w1*x1;
    %%
    % Запись полного сигнала
    xx((k-1)*length(x)+1:k*length(x)) = x;
    dd((k-1)*length(d)+1:k*length(d)) = d;
    yy((k-1)*length(y)+1:k*length(y)) = y;
end

subplot(2,1,1)
plot([xx, dd, yy]);
subplot(2,1,2)
plot(e)



% release(scope)