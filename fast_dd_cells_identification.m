clear;
close all;

pwd1 = fileparts(which('fast_dd_cells_identification.m'));

addpath([pwd1 '\general_functions']);

data_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\11_24_21_pt3';

proc_data_path = 'C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data';
% fname = 'A1_ammn1_5plt_1plm_12_27_20_OA';

fname = 'AC_rest1_mpl5';
file_num = '2';
file_date = '11_24_21';

load_dir = [data_dir '\' fname '-00' file_num];

%save_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\11_24_21_pt3';
%save_dir = 'C:\Users\ys2605\Desktop\stuff\random_save_path';
save_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data';
save_dir_movie = [save_dir '\movies'];

multiplane = 5;

mpl_tags = {'Ch2_000001', 'Ch2_000002', 'Ch2_000003', 'Ch2_000004', 'Ch2_000005'};
mpl_tags_red = {'Ch1_000001', 'Ch1_000002', 'Ch1_000003', 'Ch1_000004', 'Ch1_000005'};


%% load movie
Y = cell(multiplane,1);

for n_pl = 1:multiplane
    temp_fname = sprintf('%s_%s_pl%d.h5',fname, file_date, n_pl);
    temp_fpath = [save_dir_movie '\' temp_fname];
    if exist(temp_fpath, 'file')
        disp(['Loading ' temp_fname]);
        Y{n_pl} = h5read(temp_fpath,'/mov');
    else
        Y{n_pl} = f_collect_prairie_tiffs4(load_dir, mpl_tags{n_pl});
        f_save_mov_YS(Y{n_pl}, [save_dir_movie '\' temp_fname], '/mov');
    end
end

ave_pre = cell(multiplane,1);
for n_pl = 1:multiplane
    ave_pre{n_pl} = mean(Y{n_pl},3);
end

%% register
Y_reg = cell(multiplane,1);
dsall = cell(multiplane,1);

params_reg.image_target = [];
params_reg.num_iterations = 2;
params_reg.plot_stuff = 1;
params_reg.save_smooth = 0;
params_reg.save_reg = 0;
params_reg.save_path = save_dir_movie;
params_reg.smooth_std = [0.5 0.5 3];

fname_dsall = sprintf('%s_%s_reg_data', fname, file_date);
for n_pl = 1:multiplane
    temp_fname_reg = sprintf('%s_%s_pl%d_reg.h5',fname, file_date, n_pl);
    temp_fpath = [save_dir_movie '\' temp_fname_reg];
    if exist(temp_fpath, 'file')
        disp(['Loading ' temp_fname]);
        Y_reg{n_pl} = h5read(temp_fpath,'/mov');
        
        load_data = load([save_dir_movie '\' fname_dsall]);
        dsall{n_pl} = load_data.dsall{n_pl};
    else
        params_reg.save_fname = [fname '_' file_date '_pl' num2str(n_pl)];
        [Y_reg{n_pl}, dsall{n_pl}] = f_mpl_register2(Y{n_pl}, params_reg);

        % save reg
        f_save_mov_YS(Y_reg{n_pl}, temp_fpath, '/mov')
    end
end
clear Y;

save([save_dir_movie '\' fname_dsall], 'dsall');

ave_reg = cell(multiplane,1);
ave_reg_std = cell(multiplane,1);
for n_pl = 1:multiplane
    ave_reg{n_pl} = mean(Y_reg{n_pl},3);
    ave_reg_std{n_pl} = std(Y_reg{n_pl},[],3);
end

%% check registration
figure; 
subplot(2,1,1);
imagesc(cat(2,ave_pre{:}));
title('Before registration')
subplot(2,1,2);
imagesc(cat(2,ave_reg{:}));
title('After registration')

%%
Y_red = cell(multiplane, 1);
for n_pl = 1:multiplane
    temp_fname = sprintf('%s_%s_pl%d.h5',fname, file_date, n_pl);
    temp_fpath = [save_dir_movie '\' temp_fname];
    Y_red{n_pl} = f_collect_prairie_tiffs4(load_dir, mpl_tags_red{n_pl});
    
    temp_dsall = sum(reshape(cat(2,dsall{n_pl}{:}), [], 2, params_reg.num_iterations),3);
    Y_red{n_pl} = uint16(f_suite2p_reg_apply(Y_red{n_pl}, temp_dsall));
end

ave_red_reg = cell(multiplane,1);
for n_pl = 1:multiplane
    ave_red_reg{n_pl} = mean(Y_red{n_pl},3);
end

%% make some ave im planes
clim_prc_th = 99.8;

green_chan = cat(2,ave_reg{:});
green_chan = green_chan - min(green_chan(:));
clim_max = prctile(green_chan(:), clim_prc_th);
green_chan = green_chan/clim_max;
green_chan(green_chan>1) = 1;

red_chan = cat(2,ave_red_reg{:});
red_chan = red_chan - min(red_chan(:));
clim_max = prctile(red_chan(:), clim_prc_th);
red_chan = red_chan/clim_max;
red_chan(red_chan>1) = 1;

imall = cat(3, red_chan, green_chan, zeros(size(green_chan)));
figure; imagesc(imall/max(imall(:)));
axis equal tight
title(sprintf('%s planes 1-%d', [fname '_' file_date], multiplane), 'interpreter', 'none');

disp('saving ave images')
imwrite(imall,[save_dir_movie '\' sprintf('%s_%s_ave_im_pl_all.png', fname, file_date)])
imall_split = reshape(imall, [256, 256, 5, 3]);
for n_pl = 1:multiplane
    imwrite(squeeze(imall_split(:,:,1,:)),[save_dir_movie '\' sprintf('%s_%s_ave_im_pl%d.png', fname, file_date,n_pl)])
end

clear Y_red;
%%

fname = 'AC_ammn2_mpl5';
file_num = '2';
file_date = '11_24_21';
load_dir = [data_dir '\' fname '-00' file_num];


%%

Y = cell(multiplane,1);

for n_pl = 1:multiplane
    temp_fname = sprintf('%s_%s_pl%d.h5',fname, file_date, n_pl);
    temp_fpath = [save_dir_movie '\' temp_fname];
    if exist(temp_fpath, 'file')
        disp(['Loading ' temp_fname]);
        Y{n_pl} = h5read(temp_fpath,'/mov');
    else
        Y{n_pl} = f_collect_prairie_tiffs4(load_dir, mpl_tags{n_pl});
        f_save_mov_YS(Y{n_pl}, [save_dir_movie '\' temp_fname], '/mov');
    end
end

ave_pre_ammn = cell(multiplane,1);
for n_pl = 1:multiplane
    ave_pre_ammn{n_pl} = mean(Y{n_pl},3);
end

%%

Y_reg = cell(multiplane,1);
dsall = cell(multiplane,1);

params_reg.num_iterations = 2;
params_reg.plot_stuff = 1;
params_reg.save_smooth = 0;
params_reg.save_reg = 0;
params_reg.save_path = save_dir_movie;
params_reg.smooth_std = [0.5 0.5 3];

fname_dsall = sprintf('%s_%s_reg_data', fname, file_date);
for n_pl = 1:multiplane
    temp_fname_reg = sprintf('%s_%s_pl%d_reg.h5',fname, file_date, n_pl);
    temp_fpath = [save_dir_movie '\' temp_fname_reg];
    if exist(temp_fpath, 'file')
        disp(['Loading ' temp_fname]);
        Y_reg{n_pl} = h5read(temp_fpath,'/mov');
        
        load_data = load([save_dir_movie '\' fname_dsall]);
        dsall{n_pl} = load_data.dsall{n_pl};
    else
        params_reg.save_fname = [fname '_' file_date '_pl' num2str(n_pl)];
        params_reg.image_target = ave_reg{n_pl};
        [Y_reg{n_pl}, dsall{n_pl}] = f_mpl_register2(Y{n_pl}, params_reg);

        % save reg
        f_save_mov_YS(Y_reg{n_pl}, temp_fpath, '/mov')
    end
end
clear Y;

save([save_dir_movie '\' fname_dsall], 'dsall');

ave_reg_ammn = cell(multiplane,1);
for n_pl = 1:multiplane
    ave_reg_ammn{n_pl} = mean(Y_reg{n_pl},3);
end

%% check registration
figure; 
subplot(2,1,1);
imagesc(cat(2,ave_pre_ammn{:}));
title('Before registration')
subplot(2,1,2);
imagesc(cat(2,ave_reg_ammn{:}));
title('After registration')

%% analysis
n_pl = multiplane;
base_onset_win = [5 10];
tt = 270;

%%
proc_data = load([proc_data_path '\' fname '_processed_data']);

%%

Y2 = cat(2,Y_reg{:});

[d1, d2, T] = size(Y2);

f_save_mov_YS(Y2(:,:,10001:15000), sprintf('%s\\movie_all_mpl_reg.h5', save_dir), '/mov')

%f_save_tif_stack2_YS(Y_mpl_all, [save_path, '\' sprintf('all_pl_reg')]);

Y2 = reshape(Y2, [d1*d2, T]);
%figure; imagesc(mean(Y{n_pl},3))

stim_frame_index = proc_data.data.stim_frame_index{1};
trial_types = proc_data.data.trial_types;
MMN_orientations = proc_data.data.MMN_orientations;
stim_frame_index1 = stim_frame_index(trial_types == tt);
num_trials = numel(stim_frame_index1);

Y_sort = f_get_stim_trig_resp(Y2, stim_frame_index1, base_onset_win);
base_Y = mean(Y2,2);
base_Y2 = mean(mean(Y_sort(:,1:5,:),3),2);
Y_ave = reshape(mean(Y_sort,3)-base_Y2, d1, d2, []);%


f_save_tif_stack2_YS(Y_ave, [save_dir, '\' sprintf('mpl%d_tt%d_trial_ave',n_pl, tt)])

%%

figure; hold on;
plot(ave1)
plot(proc_data.data.file_cuts_params{1}.ave_trace)
stim_frames = zeros(numel(ave1),1);
stim_frames(stim_frame_index1) = 1;
plot(stim_frames*1600)

%%
frames_analyze = 5:15;

Y_ave_2d = reshape(Y_ave(:,:,frames_analyze), d1*d2, []);
Y_sort_2d = reshape(Y_sort(:,frames_analyze,:), d1*d2, []);

Y_in = Y_sort_2d;

method1 = 'nmf'; % 'pca', 'nmf'

if strcmpi(method1, 'pca')
    [coeff,score,latent,tsquared,explained,mu] = pca(Y_in);

elseif strcmpi(method1, 'nmf')
    [W,H] = nnmf(Y_in,4);
    coeff = H';
    score = W;
end

%%
n_pc = 1;

score_2d = reshape(score(:,n_pc),d1, d2, []);
figure; imagesc(score_2d); axis equal tight; title(['comp ' num2str(n_pc)])


coeff_sort = reshape(coeff, numel(frames_analyze), num_trials, []);
figure; plot(coeff_sort(:,:,n_pc)); title(['comp ' num2str(n_pc)])
figure; plot(mean(coeff_sort(:,:,n_pc),2)); title(['comp ' num2str(n_pc)])

%%
tr = 1:num_trials;
figure; plot3(coeff_sort(:,tr,2), coeff_sort(:,tr,3), coeff_sort(:,tr,1))
xlabel('PC2'); ylabel('PC3'); zlabel('PC1');

%%
lr_pcs = 1:3;
mov_low_rank = reshape(score(:,lr_pcs)*coeff(:,lr_pcs)',d1,d2,[]);
f_save_tif_stack2_YS(mov_low_rank, [save_dir, '\' sprintf('mpl%d_tt%d_low_rank_%s',n_pl, tt,method1)])
mov_low_rank_ave = reshape(mov_low_rank, d1, d2, numel(frames_analyze), []);
f_save_tif_stack2_YS(mov_low_rank_ave, [save_dir, '\' sprintf('mpl%d_tt%d_low_rank_ave_%s',n_pl, tt,method1)])

