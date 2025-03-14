% Deng et al., neuron, 2025 - fiber photometry analysis code
clear
close all
clc
%%
load('ff_data.mat') % load fiber photometry raw data
load('brainState.mat') % load brainstates data
ts = 0.98304; % timestamps

% data average points m
m = 1000;
fdata = ff465.data;

% down to ~1Hz
n = floor(length(fdata)/m);
temp = zeros(m,n);
temp(1:end) = fdata(1:m*n);
ff465_averaged = mean(temp);
figure
plot(ff465_averaged(5:end))
hold on

try
    ffdata = ff466.data;
    n = floor(length(ffdata)/m);
    temp = zeros(m,n);
    temp(1:end) = ffdata(1:m*n);
    ff466_averaged = mean(temp);
    %figure
    plot(ff466_averaged(5:end),'r')
    hold off
end

%% fit ff_data
ff465_averaged =  ff465_averaged';
ff466_averaged =  ff466_averaged';

p_01 = [465;100;166;166;466;100;106;106];

% process -- fitting (use NREM)
ff465_averaged(1,1) = ff465_averaged(2,1);
ff465_averaged(end,1) = ff465_averaged(end-1,1);
ff466_averaged(1,1) = ff466_averaged(2,1);
ff466_averaged(end,1) = ff466_averaged(end-1,1);

% fit - 465
ff465sc = ff465_averaged - p_01(3);
ff466sc = ff466_averaged - p_01(7);

% use NREM
ff = [];
ff(:,1) = ff465sc;
ff(:,2) = ts*(1:length(ff));
ff_sig = ff;
ff_sig = ff_brainstate(ff_sig,brainState);
idx = find(ff_sig(:,3)==0);
[fittedmodel, gof] = ffExpFit(ff_sig(idx,1),idx);
drawnow
x = 1:length(ff);
ff2= ((ff_sig(:,1) - fittedmodel(x))./fittedmodel(x));
% save
ffData.ff(1,:) = ff2;
ffData.ts = ff_sig(:,3);
% 466
ff_sig(:,1) = ff466sc;
[fittedmodel, gof] = ffExpFit(ff_sig(idx,1),idx);
drawnow
ff2= ((ff_sig(:,1) - fittedmodel(x))./fittedmodel(x));
% save
ffData.ff(2,:) = ff2;
% plot
plot(zscore(ffData.ff(1,:)))
hold on
plot(zscore(ffData.ff(2,:)),'g')
hold off
drawnow









