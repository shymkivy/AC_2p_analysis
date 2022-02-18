

clear;
close all;

fpath = 'C:\Users\ys2605\Desktop\stuff\AC_data\12_24_21b_dream\cell_traces';

traces_rest_pre = {'cell1_rest1.csv';...
                   'cell2_rest1.csv';...
                   %'cell3_rest1.csv';...
                   %'Cell4_pre.csv';...
                   'cell3_rest1.csv'};

traces_ammn_pre = {'cell1_ammn1.csv';...
                   'cell2_ammn1.csv';...
                   %'Cell3_pre.csv';...
                   %'Cell4_pre.csv';...
                   'cell3_ammn1.csv'};

traces_ammn_stim = {'cell1_ammn_stim2.csv';...
                    'cell2_ammn_stim2.csv';...
                    %'Cell3_stim.csv';...
                    %'Cell4_stim.csv';...
                    'cell3_ammn_stim2.csv'};
                 
traces_ammn_post = {'cell1_ammn_post4.csv';...
                    'cell2_ammn_post4.csv';...
                    %'Cell3_post.csv';...
                    %'Cell4_post.csv';...
                    'cell3_ammn_post4.csv'};
                 
traces_rest_post = {'cell1_rest_post3.csv';...
                    'cell2_rest_post3.csv';...
                    %'Cell3_pre.csv';...
                    %'Cell4_pre.csv';...
                    'cell3_rest_post3.csv'};
                 

cuts_stim_fname = 'AC_ammn_stim2_12_24_21b_h5cutsdata.mat';
volt_dtim_fname = 'raw_AC_ammn_stim2_12_24_21b.csv';
%%
traces_pre = cell(numel(traces_ammn_pre),1);
for n_file = 1:numel(traces_ammn_pre)
    temp_data = csvread([fpath '/' traces_ammn_pre{n_file}],1,0);
    traces_pre{n_file} = temp_data(:,2);
end
traces_pre = cat(2,traces_pre{:})';
traces_pre_d = f_smooth_dfdt3(traces_pre, 1,5);

traces_stim = cell(numel(traces_ammn_stim),1);
for n_file = 1:numel(traces_ammn_stim)
    temp_data = csvread([fpath '/' traces_ammn_stim{n_file}],1,0);
    traces_stim{n_file} = temp_data(:,2);
end
traces_stim = cat(2,traces_stim{:})';
traces_stim_d = f_smooth_dfdt3(traces_stim, 1,5);

traces_post = cell(numel(traces_ammn_post),1);
for n_file = 1:numel(traces_ammn_post)
    temp_data = csvread([fpath '/' traces_ammn_post{n_file}],1,0);
    traces_post{n_file} = temp_data(:,2);
end
traces_post = cat(2,traces_post{:})';
traces_post_d = f_smooth_dfdt3(traces_post, 1,5);

%%
cuts_stim = load([fpath '\..\' cuts_stim_fname]);

volt_data_stim = csvread([fpath '\..\' volt_dtim_fname], 1, 0);

stim_chan = 2;
figure; plot(volt_data_stim(:,stim_chan))

vid_cuts_trace = logical(cuts_stim.cuts_data{4}.vid_cuts_trace);

stim_trace_full = zeros(numel(vid_cuts_trace),1);
stim_trace_full(vid_cuts_trace) = traces_stim(1,:);

%%
%traces_stim_d = traces_post_d(:,1:6000);
%traces_post_d = traces_post_d(:,6001:end);

figure; hold on; 
subplot(3,1,1);
plot(traces_pre(1,:))
subplot(3,1,2);
plot(traces_pre(2,:))
subplot(3,1,3);
plot(traces_pre(3,:))

figure; hold on; 
plot(traces_pre_d(1,:))
plot(traces_pre_d(2,:))
plot(traces_pre_d(3,:))


figure; hold on;
plot(traces_post(1,:))
plot(traces_post(2,:))
plot(traces_post(3,:))


figure; hold on; 
plot(traces_post_d(1,:))
plot(traces_post_d(2,:))
plot(traces_post_d(3,:))


figure; hold on; 
plot(traces_stim(1,:))
plot(traces_stim(2,:))
plot(traces_stim(3,:))


%traces_pre1 = traces_pre(:,1:4000);
%traces_stim = traces_pre(:,4000:18800);
%traces_post = traces_pre(:,18801:end);

%traces_pre1_d = traces_pre_d(:,1:4000);
%traces_stim_d = traces_pre_d(:,4000:18800);
%traces_post_d = traces_pre_d(:,18801:end);
%%


SI_pre = 1-f_pdist_YS(traces_pre_d, 'cosine');
SI_stim = 1-f_pdist_YS(traces_stim_d, 'cosine');
SI_post = 1-f_pdist_YS(traces_post_d, 'cosine');

SI_ave = (SI_pre + SI_stim + SI_post)./3;

SI_all = [SI_pre;SI_stim;SI_post];
SI_all_mean_sub = [SI_pre-SI_ave;SI_stim-SI_ave;SI_post-SI_ave];
clim1 = [min(SI_all(:)) max(SI_all(:))];
clim2 = [min(SI_all_mean_sub(:)) max(SI_all_mean_sub(:))];
clim1(2) = .7;

figure; 
subplot(2,3,1); imagesc(SI_pre); caxis(clim1); title('pre')
subplot(2,3,2); imagesc(SI_stim); caxis(clim1); title('stim')
subplot(2,3,3); imagesc(SI_post); colorbar; caxis(clim1); title('post')
subplot(2,3,4); imagesc(SI_pre-SI_ave); caxis(clim2); title('pre mean sub')
subplot(2,3,5); imagesc(SI_stim-SI_ave); caxis(clim2); title('stim mean sub')
subplot(2,3,6); imagesc(SI_post-SI_ave); caxis(clim2); colorbar; title('post mean sub')


%% 
SI_pre = 1-f_pdist_YS(traces_pre_d, 'cosine');
SI_post = 1-f_pdist_YS(traces_post_d, 'cosine');

SI_ave = (SI_pre + SI_post)./2;

SI_all = [SI_pre;SI_post];
SI_all_mean_sub = [SI_pre-SI_ave;SI_post-SI_ave];
clim1 = [min(SI_all(:)) max(SI_all(:))];
clim2 = [min(SI_all_mean_sub(:)) max(SI_all_mean_sub(:))];
clim1(2) = .4;

figure; 
subplot(2,2,1); imagesc(SI_pre); caxis(clim1); title('pre')
subplot(2,2,2); imagesc(SI_post); colorbar; caxis(clim1); title('post')
subplot(2,2,3); imagesc(SI_pre-SI_ave); caxis(clim2); title('pre mean sub')
subplot(2,2,4); imagesc(SI_post-SI_ave); caxis(clim2); colorbar; title('post mean sub')



