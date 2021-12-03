clear;
close all;

pwd1 = fileparts(which('fast_dd_cells_identification.m'));

addpath([pwd1 '\general_functions']);

data_dir = 'C:\Users\shymk\Desktop\stuff\AC_data\11_24_21_pt3';

% proc_data_path = 'C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data';
% fname = 'A1_ammn1_5plt_1plm_12_27_20_OA';

fname = 'AC_rest1_mpl5';
file_num = '2';
file_date = '11_24_21';

load_dir = [data_dir '\' fname '-00' file_num];

%save_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\11_24_21_pt3';
%save_dir = 'C:\Users\ys2605\Desktop\stuff\random_save_path';
save_dir = 'C:\Users\shymk\Desktop\stuff\AC_data\caiman_data';
save_dir_movie = [save_dir '\movies'];

multiplane = 5;

mpl_tags = {'Ch2_000001', 'Ch2_000002', 'Ch2_000003', 'Ch2_000004', 'Ch2_000005'};
red_chan = {'Ch1_000001', 'Ch1_000002', 'Ch1_000003', 'Ch1_000004', 'Ch1_000005'};

%%
multiplane = 1;
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

ave_temp = cell(multiplane,1);
for n_pl = 1:multiplane
    ave_temp{n_pl} = mean(Y{n_pl},3);
end

%% register
params_reg.image_target = ave_temp;
params_reg.save_smooth = 1;
params_reg.save_path = save_dir_movie;
params_reg.save_fname = [fname '_' file_date];
Y_reg = f_mpl_register(Y, params_reg);

%%

Y = Y_reg;
dsall = cell(multiplane,1);
for n_pl = 1:multiplane       
    % smooth movie
    smooth_std = [0.5 0.5 3];
    %smooth_std = [0 0 3];
    Y_sm = f_smooth_movie(Y{n_pl}, smooth_std);

    %%
    temp_fname_sm = sprintf('%s_%s_pl%d_sm.h5',fname, file_date, n_pl);
    f_save_mov_YS(Y_sm, [save_dir_movie '\' temp_fname_sm], '/mov');

    %%
    tic;
    [~, dsall{n_pl}] = f_suite2p_register_YS(Y_sm, ave_reg{n_pl});
    toc
    clear Y_sm;

%         figure; plot(dsall{n_pl})
%         
%         temp_fname_reg = sprintf('%s_%s_pl%d_reg.h5',fname, file_date, n_pl);
%         f_save_mov_YS(Y_reg, [save_dir_movie '\' temp_fname_reg], '/mov');

end

color1 = parula(5);
figure; hold on;
for n_pl = 1:5
    plot(dsall{n_pl}(:,1), 'color', color1(n_pl,:))
end


%% apply corr to raw files
Y_reg = cell(multiplane, 1);
for n_pl = 1:multiplane
    tic;
    Y_reg{n_pl} = f_suite2p_apply_reg_YS(Y{n_pl}, dsall{n_pl});
    toc
end

ave_reg = cell(multiplane,1);
for n_pl = 1:multiplane
    ave_reg{n_pl} = mean(Y_reg{n_pl},3);
end

figure; 
subplot(2,1,1);
imagesc(cat(2,ave_reg{:}));
subplot(2,1,2);
imagesc(cat(2,ave_temp{:}));

%% save reg

for n_pl = 1:multiplane
    temp_fname_reg = sprintf('%s_%s_pl%d_reg.h5',fname, file_date, n_pl);
    f_save_mov_YS(Y_reg{n_pl}, [save_dir_movie '\' temp_fname_reg], '/mov')
end

%%

Y_cat = cat(2,Y{:,1});


f_save_mov_YS(Y_cat, sprintf('%s\\%s_%s_mpl5_cat_pre_reg.h5', save_dir_movie,fname, file_date), '/mov')

f_save_mov_YS(Y_cat_reg, sprintf('%s\\%s_%s_mpl5_cat_post_reg.h5', save_dir_movie,fname, file_date), '/mov')


ave1= squeeze(mean(mean(Y{1})));

%% registration
Y_reg = cell(multiplane,1);

save_mov = 1;

n_chan = 1;
for n_pl = 1:multiplane
    tic
    Y_reg{n_pl} = f_register_suite2p_YS(Y{n_pl,n_chan});
    toc
     
    figure;
    subplot(1,2,1); imagesc(squeeze(mean(Y{n_pl},3))); title(sprintf('pre reg mpl%d', n_pl)); axis equal tight;
    subplot(1,2,2); imagesc(squeeze(mean(Y_reg{n_pl},3))); title(sprintf('post reg mpl%d', n_pl)); axis equal tight;
     
    if save_mov
        f_save_mov_YS(Y{1,n_chan}, sprintf('%s\\%s_%s_red_mpl%d_pre_reg.h5', save_dir_movie,fname, file_date, n_pl), '/mov')
        f_save_mov_YS(Y_reg{1}, sprintf('%s\\%s_%s_red_mpl%d_post_reg.h5', save_dir_movie,fname, file_date, n_pl), '/mov')
    end
end

clear Y
%%
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

