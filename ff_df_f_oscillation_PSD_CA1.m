clear
clc
close all
%% set param
load('CAMP_tosum_n_28_0309.mat')
ff_fs = 1.0175;
toplot = [];
%%
for nn = 1:length(tosum)
    
    point1 = tosum(nn).point(1);
    point2 = tosum(nn).point(2);
    
    x1 = tosum(nn).data;
    
    data = x1;
    [b a] = butter(3,0.004*2/ff_fs,'high');
    data2 = filtfilt(b,a,data);
    
    n = length(data2);
    dataX = fft(data2);
    power_01 = 2*abs(dataX)/n;
    %         power_01 = smooth(power,20);
    hz = linspace(0,ff_fs/2,floor(n/2));
    idx =  dsearchn(hz',0.1);
    power_02 = power_01;
    power_02 = power_01./sum(power_01(1:idx));
    figure(2)
    plot(hz(1:idx),power_02(1:idx),'r','Linewidth',2)
    xlim([0 0.1])
    ylim([0 0.5])
    toplot = cat(1,toplot,power_02');
    tosum(nn).power_02 = power_02;
    toplot = cat(1,toplot,power_02');
end

%%
idx = 25;
hz = linspace(0,0.5086,122);;
figure
toplot = [];
for i = 1:length(tosum)
    power_02 = tosum(i).power_02(1:idx);
%     power_02 = smooth(power_02,3)
    toplot = cat(1,toplot,power_02');
end
meantmp = mean(toplot,1);
setmp = std(toplot,[],1)./sqrt(size(toplot,1));
shadedErrorBar(hz(1:idx),meantmp,setmp)
xlim([0 0.1])
ylim([0 0.2])