clear;
close all;

data_dir = 'G:\data\Auditory\caiman_data_missmatch\movies\';

fname = 'M1_im1_A1_ammn1_10_2_18.h5';

%%

Y = h5read([data_dir, fname], '/mov');

%%
[d1, d2, T] = size(Y);
Y2d = double(reshape(Y, [d1*d2, T]));

Ymeanim = mean(Y2d,2);
Ymeantrace = mean(Y2d);
Y2dn = Y2d - Ymeanim;

% get location of artifact
prc_thresh1 = prctile(Ymeanim, 5);
Ymask = Ymeanim<prc_thresh1;
mask_trace = mean(Y2dn(Ymask,:));

[~, idx_art]  = max(mask_trace);

%% remove components 1 by 1
num_split = 10;
update_mask = 1;
thresh_clip = 20;

Y2d_cut = Y2d(:,(idx_art-100):(idx_art+100));
[~,T2] = size(Y2d_cut);

Y2d_in = Y2d_cut;
if_plot_prc_split_mean(Y2d_in, Ymeanim, num_split);
sgtitle('before artif removal')

for n_iter = 1:10
    % normalize
    Y2d_inn = Y2d_in - mean(Y2d_in,2);
    
    % create initialize comp
    if or(n_iter == 1, update_mask)
        mask_trace_cut = mean(Y2d_inn(Ymask,:));
        mask_trace_cut(mask_trace_cut<0) = 0;
        mask_trace_cut = mask_trace_cut/norm(mask_trace_cut);
    end
    
    % compute spacial mask
    mask_spac_temp = Y2d_inn*mask_trace_cut';
    
    low_thresh = prctile(mask_spac_temp, thresh_clip);
    high_thresh = prctile(mask_spac_temp, 100 - thresh_clip);
    mask_spac_temp(mask_spac_temp<low_thresh) = low_thresh;
    mask_spac_temp(mask_spac_temp>high_thresh) = high_thresh;
    
    %[f,xi] = ksdensity(mask_spac_temp, linspace(min(mask_spac_temp(:)), max(mask_spac_temp(:)),1000), 'bandwidth', 10);
    %figure; plot(xi, f)
    
    % compute bkg;
    Ybkg_temp = mask_spac_temp*mask_trace_cut;
    
    Y2d_in = Y2d_in - Ybkg_temp;
    Y2d_in(Y2d_in<0) = 0;
    
    figure; 
    subplot(3,1,1);
    plot(mask_trace_cut); title('temporal comp')
    subplot(3,1,2:3); 
    imagesc(reshape(mask_spac_temp, d1, d2)); title('spatial comp')
    sgtitle(sprintf('inter %d', n_iter))
    if_plot_prc_split_mean(Y2d_inn, Ymeanim, num_split);
    sgtitle(sprintf('inter %d, normalized', n_iter))
end

if_plot_prc_split_mean(Y2d_in, Ymeanim, num_split);
sgtitle('after artif removal')


Ybkg = Y2d_cut - Y2d_in;

f_save_mov_YS(uint16(reshape(Y2d_cut, d1, d2, T2)), 'test_Y_cut.h5');
f_save_mov_YS(uint16(reshape(Y2d_in, d1, d2, T2)), 'test_Y_cut_dn.h5');
f_save_mov_YS(uint16(reshape(Ybkg - min(Ybkg(:)), d1, d2, T2)), 'test_Y_cut_res.h5');

%%

Y2d_out = Y2d;
Y2d_out(:,(idx_art-100):(idx_art+100)) = Y2d_out(:,(idx_art-100):(idx_art+100)) - Ybkg;
Y2d_out = uint16(round(Y2d_out));

%%
[~, fname2, ext] = fileparts(fname);

%%
f_save_mov_YS(reshape(Y2d_out, [d1, d2, T]), [data_dir, fname2, '_denoised', ext], '/mov')

%%
% 
% 
% % create initialize comp
% mask_trace_cut = mean(Y2d_cutn(Ymask,:));
% mask_trace_cutn = mask_trace_cut;
% mask_trace_cutn(mask_trace_cut<0) = 0;
% mask_trace_cutn = mask_trace_cutn/norm(mask_trace_cutn);
% 
% 
% 
% if_plot_prc_split_mean(Y2d_cutn, Ymeanim, num_split);
% 
% 
% % run nmf
% num_comp = 1;
% opt = statset('MaxIter',1000,'Display','final', 'TolFun', 1e-20, 'TolX', 1e-20);
% [W, H] = nnmf(Y2d_cutn, num_comp, 'H0', mask_trace_cutn, 'Replicates',5,...
%                    'Options',opt,...
%                    'Algorithm','mult');
% 
% %W1 = Y2d_cutn*H';
% Ybkg = W*H;
% 
% Y2d_cutdn = Y2d_cut - Ybkg;
% 
% Y2d_cutdn2 = Y2d_cutdn;
% Y2d_cutdn2(Y2d_cutdn<0) = 0;
% 
% 
% figure; subplot(4,1,1);
% plot(mean(Y2d_cutn,1)); title('artifact in normalized cut vid')
% subplot(4,1,2); 
% plot(H'); title('nmf comp of afrtifact')
% subplot(4,1,3:4);  imagesc(reshape(W(:,1), d1, d2)); title('nmf comp of afrtifact')
% 
% f_save_mov_YS(uint16(reshape(Y2d_cut, d1, d2, T2)), 'test_Y_cut.h5');
% f_save_mov_YS(uint16(reshape(Y2d_cutdn, d1, d2, T2)), 'test_Y_cut_dn.h5');
% 
% 
% if_plot_prc_split_mean(Y2d_cutdn2, Ymeanim, num_split);
% 
% 
% Y2d_cutn = Y2d_cutdn2 - mean(Y2d_cutdn2,2);
% 
% % create initialize comp
% mask_trace_cut = mean(Y2d_cutn(Ymask,:));
% mask_trace_cutn = mask_trace_cut;
% mask_trace_cutn(mask_trace_cut<0) = 0;
% mask_trace_cutn = mask_trace_cutn/norm(mask_trace_cutn);
% 
% num_split = 10;
% 
% if_plot_prc_split_mean(Y2d_cutn, Ymeanim, num_split);
% 
% 
% prc_thresh2 = prctile(Ymeanim, 5);
% Ymask2 = Ymeanim<prc_thresh2;
% Ymask_trace2 = mean(Y2d_cutn(Ymask2,:));
% Ymask_trace2(Ymask_trace2 < 0) = 0;
% Ymask_trace2 = Ymask_trace2/norm(Ymask_trace2);
% figure; plot(Ymask_trace2)
% 
% H = Ymask_trace2;
% W1 = Y2d_cutn*H';
% 
% % run nmf
% num_comp = 1;
% opt = statset('MaxIter',1000,'Display','final', 'TolFun', 1e-10, 'TolX', 1e-10);
% [W, H] = nnmf(Y2d_cutn, num_comp, 'H0', Ymask_trace2, 'Replicates',5,...
%                    'Options',opt,...
%                    'Algorithm','mult');
%            
% %W1 = Y2d_cutn*H';
% Ybkg = W1*H;
% 
% Y2d_cutdn3 = Y2d_cutdn2 - Ybkg;
% 
% figure; subplot(4,1,1);
% plot(mean(Y2d_cutn,1)); title('artifact in normalized cut vid')
% subplot(4,1,2); 
% plot(H'); title('nmf comp of afrtifact')
% subplot(4,1,3:4);  imagesc(reshape(W1(:,1), d1, d2)); title('nmf comp of afrtifact')
% 
% if_plot_prc_split_mean(Y2d_cutdn3, mean(Y2d_cutdn2,2), num_split);
% 
% 
% Y2d_cutdn4 = Y2d_cutdn3;
% Y2d_cutdn4(Y2d_cutdn3<0) = 0;
% 
% if_plot_prc_split_mean(Y2d_cutdn4, mean(Y2d_cutdn2,2), num_split);
% 
% figure; plot(std(Y2d_cutdn4))
% 
% 
% 
% f_save_mov_YS(uint16(reshape(Y2d_cut, d1, d2, T2)), 'test_Y_cut.h5');
% f_save_mov_YS(uint16(reshape(Y2d_in, d1, d2, T2)), 'test_Y_cut_dn4.h5');
% 
% 
% Y_cutdn2 = reshape(Y2d_cutdn2, d1, d2, T2);
% figure; imagesc(mean(Y_cutdn2,3))
% 
% figure; plot(squeeze(mean(mean(Y_cutdn2(230:end,:,:),1),2)))
% 
% 
% figure; plot(squeeze(Y_cutdn2(250,50:55,:))')

%%
function if_plot_prc_split_mean(data_2d, Ymeanim, num_split)

prc_range = linspace(0, 100, num_split);
prc_thresh = prctile(Ymeanim, prc_range);

figure;
for n_prc = 2:num_split
    subplot(num_split-1,1,n_prc-1)
    Ymask = and(Ymeanim>prc_thresh(n_prc-1),Ymeanim<prc_thresh(n_prc));
    mask_trace = mean(data_2d(Ymask,:));
    plot(mask_trace); axis tight;
end


end
