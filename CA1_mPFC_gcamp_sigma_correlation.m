% Deng et al., neuron, 2025 - CA1 and mpfc gcmap correlation analysis code
clear
close all
clc
%%
load('CA1_MPFC_sigma_data.mat')
% misc paramerters
ff_xcorr = [];
x_xcorr = [];
dt = 5;
cc_sum = [];


x_xcorr = [];
ff_xcorr = [];
ff = ffData.ff;  % fiber photometry data
N = floor(length(ff)/dt);
for i = 1:N-1*dt
    idx = (i-1)*dt+1:(i+1)*dt;
    temp = xcorr(ff(1,idx),ff(2,idx),'coef');
    ff_xcorr(i,:) = temp;
    x_xcorr(i) = i*dt*0.98304;
end

%% sigma phase
% sync sigma with bs
x = (1:length(ffData.bs))*0.98304;% fs = 0.98304
v = ffData.bs';
xq = (1:length(ffData.sigma))*2;% fs = 2
bs_sigma = interp1(x,v,xq)';

% band pass filter sigma
fs = 0.5;
uslowFilt = designfilt('bandpassfir','FilterOrder',500, ...
    'CutoffFrequency1',0.01,'CutoffFrequency2',0.03, ...
    'SampleRate',fs);
%fvtool(uslowFilt)
x = ffData.sigma;
xfilted = filtfilt(uslowFilt, x);
y = hilbert(xfilted);
sigphase = angle(y);

%% determine ff_xcorr phase
% sync sigma to xcorr
x = (1:length(ffData.sigma))*2;
v = sigphase;
xq = x_xcorr;
phase_xcorr = interp1(x,v,xq)';
% sigma filted
v = xfilted;
sigma_resampled = interp1(x,v,xq)';
% bs
x = (1:length(ffData.bs))*0.98304;
v = ffData.bs';
xq = x_xcorr;
bs_xcorr = interp1(x,v,xq)';


%% NREM peak or trough
dd = pi/12;
idx_peak= find(bs_xcorr == 0 & (phase_xcorr > -dd & phase_xcorr < dd));
idx_trough= find(bs_xcorr == 0 & (phase_xcorr > (pi-dd) | phase_xcorr < (-pi + dd)));
figure
plot(mean(ff_xcorr(idx_peak,:)))
hold on
plot(mean(ff_xcorr(idx_trough,:)),'r')
hold off

%% plot
subplot(411)
imagesc(ff_xcorr')
hold on
cc = ff_xcorr(:,ceil(size(ff_xcorr,2)/2));
plot(-cc*3+ceil(size(ff_xcorr,2)/2),'m')
hold off
subplot(412)
plot(x_xcorr,cc)
axis tight
subplot(413)
%plot((1:length(ffData(fN).sigma))*2,zscore(xfilted)*4,'m')
sigma_z = zscore(sigma_resampled);
plot(x_xcorr,sigma_z)
hold on
plot(x_xcorr,cc*2)
try
    plot(x_xcorr(idx_peak),sigma_z(idx_peak),'r.')
    plot(x_xcorr(idx_trough),sigma_z(idx_trough),'g.')
end
hold off
axis tight
%     temp = axis;
%     axis([temp(1) temp(2) -pi*1.3 pi*1.3])
subplot(414)
plot(x_xcorr,bs_xcorr)
axis tight

%% save
cc_sum(fN,1,:) = mean(ff_xcorr(idx_peak,:));
cc_sum(fN,2,:) = mean(ff_xcorr(idx_trough,:));


