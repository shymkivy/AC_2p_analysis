clear

data_dir = 'D:\data\AC\2p\5_21_20_im\';
fname_stim = 'A1_asynch1_5_21_20_stim_data.mat';
fname_volt = 'raw_A1_asynch1_5_21_20.csv';
fname_OA = 'F:\data\Auditory\caiman_out\OA_outputs\A1_asynch1_5_21_20_OA_cut_5gsig_results_sort.mat';

data_OA = load(fname_OA);
data_stim = load([data_dir fname_stim]);


data_volt = csvread([data_dir fname_volt], 1,0);


fp = 0.035141528;


data_stim

time_offset = data_stim.synch_pause_time(1) + data_stim.synch_pause_time(2) + fp*30;

samp1 = 195;

figure; hold on;

plot((1:(1000*samp1))/samp1, data_stim.asynch_stim(1,1:(1000*samp1)))
plot(data_stim.asynch_stim_dsp(1:1000))



figure; plot(data_volt(:,5))

sec_ave = .5;

figure;
%spectrogram(data_stim.asynch_stim(1,1:(2000*samp1)))

spectrogram(data_stim.asynch_stim(1:round(sec_ave/data_stim.sig_dt)),1000,100,1000, 1/data_stim.sig_dt, 'yaxis');

15194
625812

figure; plot(data_OA.est.S(3,:))