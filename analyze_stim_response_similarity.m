

clear;
close all;

fpath = 'C:\Users\ys2605\Desktop\stuff\AC_data\11_24_21_pt3\cell_traces';

traces_names_pre = {'Cell1_pre.csv';...
                    'Cell2_pre.csv';...
                    'Cell3_pre.csv';...
                    'Cell4_pre.csv';...
                    'Cell5_pre.csv'};

traces_names_stim = {'Cell1_stim.csv';...
                     'Cell2_stim.csv';...
                     'Cell3_stim.csv';...
                     'Cell4_stim.csv';...
                     'Cell5_stim.csv'};
                 
traces_names_post = {'Cell1_post.csv';...
                     'Cell2_post.csv';...
                     'Cell3_post.csv';...
                     'Cell4_post.csv';...
                     'Cell5_post.csv'};
                 

cuts_stim_fname = 'AC_stim2_mpl5_11_24_21_h5cutsdata.mat';
volt_dtim_fname = 'raw_AC_stim2_11_24_21_mpl5.csv';
%%
traces_pre = cell(numel(traces_names_pre),1);
for n_file = 1:numel(traces_names_pre)
    temp_data = csvread([fpath '/' traces_names_pre{n_file}],1,0);
    traces_pre{n_file} = temp_data(:,2);
end
traces_pre = cat(2,traces_pre{:})';
traces_pre_d = f_smooth_dfdt3(traces_pre, 1,5);

traces_stim = cell(numel(traces_names_stim),1);
for n_file = 1:numel(traces_names_stim)
    temp_data = csvread([fpath '/' traces_names_stim{n_file}],1,0);
    traces_stim{n_file} = temp_data(:,2);
end
traces_stim = cat(2,traces_stim{:})';
traces_stim_d = f_smooth_dfdt3(traces_stim, 1,5);

traces_post = cell(numel(traces_names_post),1);
for n_file = 1:numel(traces_names_post)
    temp_data = csvread([fpath '/' traces_names_post{n_file}],1,0);
    traces_post{n_file} = temp_data(:,2);
end
traces_post = cat(2,traces_post{:})';
traces_post_d = f_smooth_dfdt3(traces_post, 1,5);

%%
cuts_stim = load([fpath '\' cuts_stim_fname]);



volt_data_stim = csvread([fpath '\' volt_dtim_fname], 1, 0);

stim_chan = 6;
figure; plot(volt_data_stim(:,stim_chan))

vid_cuts_trace = cuts_stim.cuts_data{3}.vid_cuts_trace;

stim_trace_full = zeros(numel(vid_cuts_trace),1);
stim_trace_full(vid_cuts_trace) = traces_stim(1,:);

%%
%traces_stim_d = traces_post_d(:,1:6000);
%traces_post_d = traces_post_d(:,6001:end);

figure; hold on; 
subplot(2,1,1);
plot(traces_pre(1,:))
subplot(2,1,2);
plot(traces_pre(2,:))

figure; hold on; 
plot(traces_pre_d(1,:))
plot(traces_pre_d(2,:))

figure; hold on;
plot(traces_post(1,:))
plot(traces_post(2,:))

figure; hold on; 
plot(traces_post_d(1,:))
plot(traces_post_d(2,:))

figure; hold on; 
plot(traces_stim(1,:))
plot(traces_stim(2,:))


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



