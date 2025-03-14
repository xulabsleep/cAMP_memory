% Deng et al., neuron, 2025 - CAMP phase of ripple and spindle
clear
close all
clc
%% load data 
load('ripple_spindle_dataNex.mat');  % ripple and spindle data

% process event
contvars_list = {};
for i = 1:length(dataNex.contvars)
    contvars_list{i} = dataNex.contvars{i,1}.name;
end

temp = strfind(contvars_list,'fff46');
temp_idx = find(not(cellfun('isempty',temp)));
fff465_01 = dataNex.contvars{temp_idx,1}.data;
%     fff465_01 = smooth(fff465_01,10,'sgolay');
ff_fs = dataNex.contvars{temp_idx, 1}.ADFrequency;
ts = 1./ff_fs;
fff465_01(:,2) = ts:ts:ts*length(fff465_01);
fff465_ts = timeseries(fff465_01(:,1),fff465_01(:,2));
% brainstates
temp = strfind(contvars_list,'brainState_01');
temp_idx = find(not(cellfun('isempty',temp)));
bs = dataNex.contvars{temp_idx,1}.data;
ff465bs = ff_brainstate(fff465_01,bs);
bs_ff = ff465bs(:,3);

events_list = {};
for i = 1:length(dataNex.events)
    events_list{i} = dataNex.events{i,1}.name;
end

temp = strfind(events_list,'r_peaks_nrem');
idx = find(not(cellfun('isempty',temp)));
r_peaks = dataNex.events{idx(1),1}.timestamps;
r_idx = round(r_peaks*ff_fs);

temp = strfind(events_list,'spindles_eeg_nrem');
idx = find(not(cellfun('isempty',temp)));
spindles = dataNex.events{idx(1),1}.timestamps;
spindle_idx = ceil(spindles*ff_fs);

[~,idx_peaks,~] = findpeaks(zscore(fff465_01(:,1)),...
    'MinPeakProminence',0.5,'Annotate','extents');
[~,idx_trough,~] = findpeaks(-zscore(fff465_01(:,1)),...
    'MinPeakProminence',0.5,'Annotate','extents');

peak_tmp = zeros(length(fff465_01),1);
peak_tmp(idx_peaks) = 1;

idx_peak_tmp = peak_tmp & (bs_ff==0);
idx_peaks = find(idx_peak_tmp);

figure
hold on
plot(fff465_01(:,1))
scatter(spindle_idx,fff465_ts.data(spindle_idx),15,'g','filled')
scatter(r_idx,fff465_ts.data(r_idx),15,'r','filled')

% ripple
phase_r = [];
for nn = 1:length(idx_peaks)
    idx1 = idx_peaks(nn);
    tmp = dsearchn(idx_trough,idx1);
    
    try
        if idx_trough(tmp) < idx1
            idx0 = idx_trough(tmp);
            idx2 = idx_trough(tmp+1);
        else
            if tmp == 1
                continue
            end
            idx2 = idx_trough(tmp);
            idx0 = idx_trough(tmp-1);
        end
        
        x = fff465_01(idx0:idx2,1);
        
        r_idx_tmp = r_idx(find(r_idx>idx0 & r_idx<idx2));
        r_idx_tmp = r_idx_tmp-idx0+1;
        
        H1 = fff465_01(idx1,1)-fff465_01(idx0,1);
        H2 = fff465_01(idx1,1)-fff465_01(idx2,1);
        
        phase_tmp = [];
        for xx = 1:length(x)
            if (xx-1)<(idx1-idx0) | (xx-1) == (idx1-idx0)
                phase_tmp(xx) = ((x(xx)-fff465_01(idx0,1))./H1)*pi-pi;
            else
                phase_tmp(xx) = ((H2-(x(xx)-fff465_01(idx2,1)))./H2)*pi;
            end
        end
        
        phase_r = cat(2,phase_r,phase_tmp(r_idx_tmp));
    end
end


% spindle
phase_spindle = [];
for nn = 1:length(idx_peaks)
    idx1 = idx_peaks(nn);
    tmp = dsearchn(idx_trough,idx1);
    
    try
        if idx_trough(tmp) < idx1
            idx0 = idx_trough(tmp);
            idx2 = idx_trough(tmp+1);
        else
            if tmp == 1
                continue
            end
            idx2 = idx_trough(tmp);
            idx0 = idx_trough(tmp-1);
        end
        
        x = fff465_01(idx0:idx2,1);
        
        r_idx_tmp = spindle_idx(find(spindle_idx>idx0 & spindle_idx<idx2));
        r_idx_tmp = r_idx_tmp-idx0+1;
        
        H1 = fff465_01(idx1,1)-fff465_01(idx0,1);
        H2 = fff465_01(idx1,1)-fff465_01(idx2,1);
        
        phase_tmp = [];
        for xx = 1:length(x)
            if (xx-1)<(idx1-idx0) | (xx-1) == (idx1-idx0)
                phase_tmp(xx) = ((x(xx)-fff465_01(idx0,1))./H1)*pi-pi;
            else
                phase_tmp(xx) = ((H2-(x(xx)-fff465_01(idx2,1)))./H2)*pi;
            end
        end
        phase_spindle = cat(2,phase_spindle,phase_tmp(r_idx_tmp));
    end
end

figure
h = polarhistogram(phase_r,25);
hold on
h = polarhistogram(phase_spindle,25);
ax = gca;
ax.ThetaTick = [0:30:360];
ax.ThetaTickLabel = {'0','30','60','90','120','150','180/-180',...
    '-150','-120','-90','-60','-30','0'}
legend({'ripple','spindle'})

