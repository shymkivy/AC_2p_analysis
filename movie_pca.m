clear;
close all;


pwd1 = fileparts(which('fast_dd_cells_identification.m'));

addpath([pwd1 '\general_functions']);

mov_load_path = 'E:\data\abl\3_27_21_abl_prairie1\fov1_before-001';

% proc_data_path = 'C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data';
% fname = 'A1_ammn1_5plt_1plm_12_27_20_OA';
%proc_data_path = 'E:\data\AC\AC_data_OA_3_16_20\';
%fname = 'DF_ammn1_10_16_18_OA';

save_path = 'C:\Users\ys2605\Desktop\stuff\random_save_path';

% multiplane = 5;
% mpl_tags = {'Ch2_000001', 'Ch2_000002', 'Ch2_000003', 'Ch2_000004', 'Ch2_000005'};

multiplane = 1;
mpl_tags = {'Ch2'};
%% load movie
Y = cell(multiplane,1);
for n_pl = 1:multiplane
    Y{n_pl} = f_collect_prairie_tiffs4(mov_load_path, mpl_tags{n_pl});
end

%proc_data = load([proc_data_path '\' fname '_processed_data']);


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
Y2 = cat(2,Y_reg{:});

[d1, d2, T] = size(Y2);

Y2 = reshape(Y2, [d1*d2, T]);

method1 = 'nmf'; % 'pca', 'nmf'

if strcmpi(method1, 'pca')
    [coeff,score,latent,tsquared,explained,mu] = pca(Y2);

elseif strcmpi(method1, 'nmf')
    [W,H] = nnmf(Y2,10);
    coeff = H';
    score = W;
end


%%


n_pc = 2;

score2 = score(:,n_pc);
z1 = std(score2);
idx1 = score2<3*z1;
score2(idx1) = 0;
score_2d = reshape(score2,d1, d2, []);
figure; imagesc(score_2d); axis equal tight; title(['comp ' num2str(n_pc)])


coeff_sort = reshape(coeff, numel(frames_analyze), num_trials, []);
figure; plot(coeff_sort(:,:,n_pc)); title(['comp ' num2str(n_pc)])
figure; plot(mean(coeff_sort(:,:,n_pc),2)); title(['comp ' num2str(n_pc)])



%%

lr_pcs = 1:10;
mov_low_rank = reshape(score(:,lr_pcs)*coeff(:,lr_pcs)',d1,d2,[]);
f_save_tif_stack2_YS(mov_low_rank, [save_path, '\' sprintf('low_rank')])





