clear;
close all;


pwd1 = fileparts(which('fast_dd_cells_identification.m'));

addpath([pwd1 '\general_functions']);

mov_load_path = 'C:\Users\ys2605\Desktop\stuff\AC_data\3_27_21_mpl\A2_ammn3_abl-001';

% proc_data_path = 'C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data';
% fname = 'A1_ammn1_5plt_1plm_12_27_20_OA';
proc_data_path = 'C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data\';
fname = 'A2_ammn3_abl1_3_27_21';

save_path = 'C:\Users\ys2605\Desktop\stuff\random_save_path';

multiplane = 5;
mpl_tags = {'Ch2_000001', 'Ch2_000002', 'Ch2_000003', 'Ch2_000004', 'Ch2_000005'};

% multiplane = 1;
% mpl_tags = {'Ch2'};
%% load movie
Y = cell(multiplane,1);
for n_pl = 1:multiplane
    Y{n_pl} = f_collect_prairie_tiffs4(mov_load_path, mpl_tags{n_pl});
end

proc_data = load([proc_data_path '\' fname '_processed_data']);


%% registration
Y_reg = cell(multiplane,1);

save_mov = 1;

for n_pl = 1:multiplane
    tic
    Y_reg{n_pl} = f_register_suite2p_YS(Y{n_pl});
    toc
     
    figure;
    subplot(1,2,1); imagesc(squeeze(mean(Y{n_pl},3))); title(sprintf('pre reg mpl%d', n_pl)); axis equal tight;
    subplot(1,2,2); imagesc(squeeze(mean(Y_reg{n_pl},3))); title(sprintf('post reg mpl%d', n_pl)); axis equal tight;
     
    if save_mov
        f_save_mov_YS(Y{1}, sprintf('%s\\movie_mpl%d_pre_reg.h5', save_path, n_pl), '/mov')
        f_save_mov_YS(Y_reg{1}, sprintf('%s\\movie_mpl%d_post_reg.h5', save_path, n_pl), '/mov')
    end
end

clear Y

%%

n_pl = 0;
base_onset_win = [5 10];
tt = 270;

%%

Y2 = cat(2,Y_reg{:});

[d1, d2, T] = size(Y2);

f_save_mov_YS(Y2(:,:,10001:15000), sprintf('%s\\movie_all_mpl_reg.h5', save_path), '/mov')

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


f_save_tif_stack2_YS(Y_ave, [save_path, '\' sprintf('mpl%d_tt%d',n_pl, tt)])



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


n_pc = 1;

score_2d = reshape(score(:,n_pc),d1, d2, []);
figure; imagesc(score_2d); axis equal tight; title(['comp ' num2str(n_pc)])


coeff_sort = reshape(coeff, numel(frames_analyze), num_trials, []);
figure; plot(coeff_sort(:,:,n_pc)); title(['comp ' num2str(n_pc)])
figure; plot(mean(coeff_sort(:,:,n_pc),2)); title(['comp ' num2str(n_pc)])


tr = 1:num_trials;
figure; plot3(coeff_sort(:,tr,2), coeff_sort(:,tr,3), coeff_sort(:,tr,1))
xlabel('PC2'); ylabel('PC3'); zlabel('PC1');


lr_pcs = 1:3;
mov_low_rank = reshape(score(:,lr_pcs)*coeff(:,lr_pcs)',d1,d2,[]);
f_save_tif_stack2_YS(mov_low_rank, [save_path, '\' sprintf('mpl%d_tt%d_low_rank',n_pl, tt)])
mov_low_rank_ave = reshape(mov_low_rank, d1, d2, numel(frames_analyze), []);
f_save_tif_stack2_YS(mov_low_rank_ave, [save_path, '\' sprintf('mpl%d_tt%d_low_rank_ave',n_pl, tt)])

