clear;
close all;

pwd1 = fileparts(which('fast_dd_cells_identification.m'));
addpath([pwd1 '\general_functions']);
proc_data_path = 'F:\AC_data\caiman_data_dream3';
save_dir_movie = [proc_data_path '\movies'];


%%

data_dir = 'F:\AC_data\M166_6_20_22_pt2_dream';
fname = 'M166_im2_AC_ammn2_6_20_22_pt2_mpl5';
 
multiplane = 5;


% 
% fname = 'AC_ammn1';
% file_num = '';
% file_date = '12_24_21a';
% 
% load_dir = [data_dir '\' fname '-00' file_num];
% 
% %save_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\11_24_21_pt3';
% %save_dir = 'C:\Users\ys2605\Desktop\stuff\random_save_path';
% save_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data';

% multiplane = 5;
% 
% mpl_tags = {'Ch2_000001', 'Ch2_000002', 'Ch2_000003', 'Ch2_000004', 'Ch2_000005'};
% mpl_tags_red = {'Ch1_000001', 'Ch1_000002', 'Ch1_000003', 'Ch1_000004', 'Ch1_000005'};


%% load movie
Y = cell(multiplane,1);

for n_pl = 1:multiplane
    temp_fname = sprintf('%s_pl%d.h5',fname, n_pl);
    temp_fpath = [save_dir_movie '\' temp_fname];
    Y{n_pl} = h5read(temp_fpath,'/mov');
end

ave_proj = cell(multiplane,1);
for n_pl = 1:multiplane
    ave_proj{n_pl} = mean(Y{n_pl},3);
end


ave_proj2 = cat(2,ave_proj{:});

for n_pl = 1:multiplane
    temp_im = ave_proj{n_pl};
    temp_im = temp_im - min(temp_im(:));
    temp_im = temp_im/max(temp_im(:));
    imwrite(temp_im,sprintf('%s\\%s_pl%d_ave_im.png', data_dir, fname, n_pl));
end

std_proj = cell(multiplane,1);
for n_pl = 1:multiplane
    std_proj{n_pl} = std(double(Y{n_pl}),[],3);
end
std_proj2 = cat(2,std_proj{:});



%% check registration
f1 = figure; 
imagesc(ave_proj2);
axis equal tight;
title('ave proj')
f1.Position(3) = 2000;
saveas(f1, [data_dir '\ave_proj.fig'])
saveas(f1, [data_dir '\ave_proj.png'])

f1 = figure; 
imagesc(std_proj2);
axis equal tight;
title('std proj')
f1.Position(3) = 2000;
saveas(f1, [data_dir '\std_proj.fig'])
saveas(f1, [data_dir '\std_proj.png'])

%%
% Y_red = cell(multiplane, 1);
% for n_pl = 1:multiplane
%     temp_fname = sprintf('%s_%s_pl%d.h5',fname, file_date, n_pl);
%     temp_fpath = [save_dir_movie '\' temp_fname];
%     Y_red{n_pl} = f_collect_prairie_tiffs4(load_dir, mpl_tags_red{n_pl});
%     
%     temp_dsall = sum(reshape(cat(2,dsall{n_pl}{:}), [], 2, params_reg.num_iterations),3);
%     Y_red{n_pl} = uint16(f_suite2p_reg_apply(Y_red{n_pl}, temp_dsall));
% end
% 
% ave_red_reg = cell(multiplane,1);
% for n_pl = 1:multiplane
%     ave_red_reg{n_pl} = mean(Y_red{n_pl},3);
% end

%% make some ave im planes
% clim_prc_th = 99.8;
% 
% green_chan = cat(2,ave_reg{:});
% green_chan = green_chan - min(green_chan(:));
% clim_max = prctile(green_chan(:), clim_prc_th);
% green_chan = green_chan/clim_max;
% green_chan(green_chan>1) = 1;
% 
% red_chan = cat(2,ave_red_reg{:});
% red_chan = red_chan - min(red_chan(:));
% clim_max = prctile(red_chan(:), clim_prc_th);
% red_chan = red_chan/clim_max;
% red_chan(red_chan>1) = 1;
% 
% imall = cat(3, red_chan, green_chan, zeros(size(green_chan)));
% figure; imagesc(imall/max(imall(:)));
% axis equal tight
% title(sprintf('%s planes 1-%d', [fname '_' file_date], multiplane), 'interpreter', 'none');
% 
% disp('saving ave images')
% imwrite(imall,[save_dir_movie '\' sprintf('%s_%s_ave_im_pl_all.png', fname, file_date)])
% imall_split = reshape(imall, [256, 256, 5, 3]);
% for n_pl = 1:multiplane
%     imwrite(squeeze(imall_split(:,:,1,:)),[save_dir_movie '\' sprintf('%s_%s_ave_im_pl%d.png', fname, file_date,n_pl)])
% end
% 
% clear Y_red;


%%
proc_data = load([proc_data_path '\preprocessing\' fname(1:end-5) '_processed_data']);
cuts_data = load([proc_data_path '\preprocessing\' fname(1:end-5) '_h5cutsdata']);

%% insert cuts

[d1, d2, ~] = size(Y{1});

for n_pl = 1:multiplane
    Y_full = zeros(d1, d2, numel(cuts_data.cuts_data{n_pl}.ave_trace), 'uint16');
    Y_full(:,:, logical(cuts_data.cuts_data{n_pl}.vid_cuts_trace)) = Y{n_pl};
    Y{n_pl} = Y_full;
end

%% analysis
n_pl = multiplane;

base_onset_win = [3 12];

tt_pairs = {170, 270;...
            3, 5};

%%

Y2 = cat(2,Y{:});

[d1, d2, T] = size(Y2);

%f_save_mov_YS(Y2(:,:,10001:15000), sprintf('%s\\movie_all_mpl_reg.h5', save_dir), '/mov')

%f_save_tif_stack2_YS(Y_mpl_all, [save_path, '\' sprintf('all_pl_reg')]);

Y2 = reshape(Y2, [d1*d2, T]);
%figure; imagesc(mean(Y{n_pl},3))

stim_frame_index = proc_data.data.stim_times_frame{1,1};%proc_data.data.stim_frame_index{1};
trial_types = proc_data.data.trial_types;
MMN_orientations = proc_data.data.MMN_orientations;

sig_z_thresh =2.5;

%%   


for n_pairs = 1:size(tt_pairs,2)
    ave_resp_on_all = cell(2,1);
    ave_resp_off_all = cell(2,1);
    for n_tt = 1:2
        tt1 = tt_pairs{n_pairs,n_tt};
        stim_frame_index1 = stim_frame_index(trial_types == tt1);
        num_trials = numel(stim_frame_index1);

        Y_sort = f_get_stim_trig_resp(Y2, stim_frame_index1, base_onset_win);
        %base_Y = mean(Y2,2);
        base_Y2 = mean(mean(Y_sort(:,1:base_onset_win(1),:),3),2);
        Y_ave = reshape(mean(Y_sort,3)-base_Y2, d1, d2, []);%

        f_save_tif_stack2_YS(Y_ave, [data_dir, '\' sprintf('mpl%d_tt%d_trial_ave',n_pl, tt1)])

        ave_resp_on = mean(Y_ave(:,:,5:9),3);
        ave_resp2_on = ave_resp_on - mean(ave_resp_on(:));
        ave_resp2_on = ave_resp2_on/std(ave_resp2_on(:));
        ave_resp_on_all{n_tt} = ave_resp2_on>sig_z_thresh;

        ave_resp_off = mean(Y_ave(:,:,9:end),3);
        ave_resp2_off = ave_resp_off - mean(ave_resp_off(:));
        ave_resp2_off = ave_resp2_off/std(ave_resp2_off(:));
        ave_resp_off_all{n_tt} = ave_resp2_off>sig_z_thresh;
    end

    ave_proj3 = ave_proj2 - min(ave_proj2(:));
    ave_proj3 = ave_proj3/max(ave_proj3(:));
    ave_proj4 = cat(3, ave_resp_on_all{1}, ave_proj3, ave_resp_on_all{2});
    
    f2 = figure;
    imagesc(ave_proj4);
    axis equal tight;
    title(sprintf('onset resp red=%d, blue=%d',tt_pairs{n_pairs,1}, tt_pairs{n_pairs,2}));
    f2.Position(3) = 2000;
    saveas(f2, sprintf('%s\\onset_resp_ttR%d_ttB%d.fig', data_dir, tt_pairs{n_pairs,1}, tt_pairs{n_pairs,2}))
    saveas(f2, sprintf('%s\\onset_resp_ttR%d_ttB%d.png', data_dir, tt_pairs{n_pairs,1}, tt_pairs{n_pairs,2}))
    
    ave_proj4 = cat(3, ave_resp_off_all{1}, ave_proj3, ave_resp_off_all{2});
    f3 = figure;
    imagesc(ave_proj4);
    axis equal tight;
    title(sprintf('offset resp red=%d, blue=%d',tt_pairs{n_pairs,1}, tt_pairs{n_pairs,2}));
    f3.Position(3) = 2000;
    saveas(f3, sprintf('%s\\offset_resp_ttR%d_ttB%d.fig', data_dir, tt_pairs{n_pairs,1}, tt_pairs{n_pairs,2}))
    saveas(f3, sprintf('%s\\offset_resp_ttR%d_ttB%d.png', data_dir, tt_pairs{n_pairs,1}, tt_pairs{n_pairs,2}))
    
end


%%
% 
% figure; hold on;
% plot(proc_data.data.file_cuts_params{1}.ave_trace)
% stim_frames = zeros(numel(proc_data.data.file_cuts_params{1}.ave_trace),1);
% stim_frames(stim_frame_index1) = 1;
% plot(stim_frames*1600)
% 
% %%
% frames_analyze = 3:15;
% 
% Y_ave_2d = reshape(Y_ave(:,:,frames_analyze), d1*d2, []);
% Y_sort_2d = reshape(Y_sort(:,frames_analyze,:), d1*d2, []);
% 
% Y_in = Y_sort_2d;
% 
% method1 = 'nmf'; % 'pca', 'nmf'
% 
% if strcmpi(method1, 'pca')
%     [coeff,score,latent,tsquared,explained,mu] = pca(Y_in);
% 
% elseif strcmpi(method1, 'nmf')
%     [W,H] = nnmf(Y_in,4);
%     coeff = H';
%     score = W;
% end
% 
% %%
% n_pc = 3;
% 
% score_2d = reshape(score(:,n_pc),d1, d2, []);
% figure; imagesc(score_2d); axis equal tight; title(['comp ' num2str(n_pc)])
% 
% 
% coeff_sort = reshape(coeff, numel(frames_analyze), num_trials, []);
% figure; plot(coeff_sort(:,:,n_pc)); title(['comp all tr ' num2str(n_pc)])
% figure; plot(mean(coeff_sort(:,:,n_pc),2)); title(['comp trial ave' num2str(n_pc)])
% 
% %%
% tr = 1:num_trials;
% figure; plot3(coeff_sort(:,tr,2), coeff_sort(:,tr,3), coeff_sort(:,tr,1))
% xlabel('PC2'); ylabel('PC3'); zlabel('PC1');
% 
% %%
% lr_pcs = 1:3;
% mov_low_rank = reshape(score(:,lr_pcs)*coeff(:,lr_pcs)',d1,d2,[]);
% f_save_tif_stack2_YS(mov_low_rank, [data_dir, '\' sprintf('mpl%d_tt%d_low_rank_%s',n_pl, tt,method1)])
% mov_low_rank_ave = reshape(mov_low_rank, d1, d2, numel(frames_analyze), []);
% f_save_tif_stack2_YS(mov_low_rank_ave, [data_dir, '\' sprintf('mpl%d_tt%d_low_rank_ave_%s',n_pl, tt,method1)])
% 
