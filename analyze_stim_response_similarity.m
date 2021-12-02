

clear;
%close all;

fpath = 'C:\Users\ys2605\Desktop\stuff\AC_data\11_10_21\rest_stim_rest_fov2_4cells-001';

traces_names_pre = {'cell1ll_trace.csv';...
                    'cell2ul_trace.csv';...
                    'cell3ur_trace.csv';...
                    'cell4lr_trace.csv'};

traces_names_post = {'cell1ll_trace.csv';...
                'cell2ul_trace.csv';...
                'cell3ur_trace.csv';...
                'cell4lr_trace.csv'};

%%
traces_pre = cell(numel(traces_names_pre),1);
for n_file = 1:numel(traces_names_pre)
    temp_data = csvread([fpath '/' traces_names_pre{n_file}],1,0);
    traces_pre{n_file} = temp_data(:,2);
end

traces_pre = cat(2,traces_pre{:})';
traces_pre_d = f_smooth_dfdt3(traces_pre, 1,5);

traces_post = cell(numel(traces_names_post),1);
for n_file = 1:numel(traces_names_post)
    temp_data = csvread([fpath '/' traces_names_post{n_file}],1,0);
    traces_post{n_file} = temp_data(:,2);
end
traces_post = cat(2,traces_post{:})';
traces_post_d = f_smooth_dfdt3(traces_post, 1,5);

traces_stim_d = traces_post_d(:,1:6000);
traces_post_d = traces_post_d(:,6001:end);

figure; hold on; 
plot(traces_pre(1,:))
plot(traces_pre(2,:))
plot(traces_pre(3,:))

figure; hold on; 
plot(traces_pre_d(1,:))
plot(traces_pre_d(2,:))
plot(traces_pre_d(3,:))

traces_pre1 = traces_pre(:,1:4000);
traces_stim = traces_pre(:,4000:18800);
traces_post = traces_pre(:,18801:end);


traces_pre1_d = traces_pre_d(:,1:4000);
traces_stim_d = traces_pre_d(:,4000:18800);
traces_post_d = traces_pre_d(:,18801:end);
%%

figure; plot(traces_stim_d(1,:))


SI_pre = 1-f_pdist_YS(traces_pre1_d, 'cosine');
SI_stim = 1-f_pdist_YS(traces_stim_d, 'cosine');
SI_post = 1-f_pdist_YS(traces_post_d, 'cosine');

SI_ave = (SI_pre + SI_stim + SI_post)./3;

SI_all = [SI_pre;SI_stim;SI_post];
SI_all_mean_sub = [SI_pre-SI_ave;SI_stim-SI_ave;SI_post-SI_ave];
clim1 = [min(SI_all(:)) max(SI_all(:))];
clim2 = [min(SI_all_mean_sub(:)) max(SI_all_mean_sub(:))];
clim1(2) = .4;

figure; 
subplot(2,3,1); imagesc(SI_pre); caxis(clim1); title('pre')
subplot(2,3,2); imagesc(SI_stim); caxis(clim1); title('stim')
subplot(2,3,3); imagesc(SI_post); colorbar; caxis(clim1); title('post')
subplot(2,3,4); imagesc(SI_pre-SI_ave); caxis(clim2); title('pre mean sub')
subplot(2,3,5); imagesc(SI_stim-SI_ave); caxis(clim2); title('stim mean sub')
subplot(2,3,6); imagesc(SI_post-SI_ave); caxis(clim2); colorbar; title('post mean sub')



