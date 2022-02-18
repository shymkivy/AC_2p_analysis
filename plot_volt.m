
fpath = 'D:\data\AC\1_2_22b_dream\AC_ammn_stim-003\';
fname = 'AC_ammn_stim-003_Cycle00001_VoltageRecording_001.csv';

data = csvread([fpath fname], 1);


figure; hold on;
plot(data(:,1)/1000, data(:,2))
plot(data(:,1)/1000, data(:,6))
plot(data(:,1)/1000, data(:,7))