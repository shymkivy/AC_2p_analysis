function f_mpl_plot_spatial_tunning(cdata, ops, cond_name, n_dset)
  
%[num_cells, n_dset] = max(cdata.num_cells);
num_cells = cdata.num_cells(n_dset);
A = cdata.OA_data{n_dset}.est.A(:,cdata.OA_data{n_dset}.proc.comp_accepted);
A = A/max(A(:));
%A = reshape(A,256,256,num_cells);

figure; imagesc(reshape(mean(A,2),256,256));
axis equal tight;
title(sprintf('%s %s %s, %d cells', cond_name,ops.paradigm_type, ops.file_names.(cond_name){n_dset},num_cells))


%% freq tuning plots
jet_map = jet(10);

rc_freq_onset = cdata.peak_tuned_trials_onset{n_dset};
rc_freq_onset_mag = cdata.tuning_freq{n_dset}.onset_tunning_mag;

if_plot_tuning_map(A,rc_freq_onset,rc_freq_onset_mag,jet_map);
title([cond_name ' Onset tuned cells n = ' num2str(sum(logical(sum(rc_freq_onset,2))))]);

rc_freq_offset = cdata.resp_cells_freq_offset{n_dset};
rc_freq_offset_mag = cdata.tuning_freq{n_dset}.offset_tunning_mag;

if_plot_tuning_map(A,rc_freq_offset,rc_freq_offset_mag,jet_map);
title([cond_name ' Offset tuned cells n = ' num2str(sum(logical(sum(rc_freq_offset,2))))]);

rc_both = logical(rc_freq_onset + rc_freq_offset);
rc_mag_both = (rc_freq_onset.*rc_freq_onset_mag + rc_freq_offset.*rc_freq_offset_mag);

if_plot_tuning_map(A,rc_both,rc_mag_both,jet_map);
title([cond_name ' Onset and Offset tuned cells n = ' num2str(sum(logical(sum(rc_both,2))))]);
%figure; imagesc(1,ops.control_carrier_freq/1000, reshape(jet(100),100,1,3));

%% deviance plot
dd_map = [.8 .8 .8; .2 .6 1; 1 .2 .2];
dd_map1 = [dd_map; dd_map];


rc_mmn_onset = cdata.resp_cells_mmn_onset{n_dset};
rc_mmn_onset_mag = cdata.tuning_all{n_dset}.onset_tunning_mag(:,cdata.ctx_mmn{n_dset});

if_plot_tuning_map(A,rc_mmn_onset,rc_mmn_onset_mag,dd_map1);
title([cond_name ' Onset context tuned cells n = ' num2str(sum(logical(sum(rc_mmn_onset,2))))]);

rc_mmn_offset = cdata.resp_cells_mmn_offset{n_dset};
rc_mmn_offset_mag = cdata.tuning_all{n_dset}.offset_tunning_mag(:,cdata.ctx_mmn{n_dset});

if_plot_tuning_map(A,rc_mmn_offset,rc_mmn_offset_mag,dd_map1);
title([cond_name ' Offset context tuned cells n = ' num2str(sum(logical(sum(rc_mmn_offset,2))))]);

rc_mmn_both = logical(rc_mmn_onset + rc_mmn_offset);
rc_mmn_both_mag = rc_mmn_onset.*rc_mmn_onset_mag + rc_mmn_offset.*rc_mmn_offset_mag;
if_plot_tuning_map(A,rc_mmn_both,rc_mmn_both_mag,dd_map1);
title([cond_name ' On and offset context tuned cells n = ' num2str(sum(logical(sum(rc_mmn_both,2))))]);
%figure; imagesc(1,ops.control_carrier_freq/1000, reshape(dd_map1,6,1,3));

end

function if_plot_tuning_map(A,resp_cells_freq,resp_cells_freq_mag, colormap1)
A = A/(max(A(:)));
[pix, num_cells] = size(A);

im_loc = zeros(pix,3);
for n_cell = 1:num_cells
    if sum(resp_cells_freq(n_cell,:))>0
        [~, freq1] = max(resp_cells_freq_mag(n_cell,:));
        temp = A(:,n_cell)>0;
        for cc = 1:3
            im_loc(temp,cc) = A(temp,n_cell)*colormap1(freq1,cc);
        end
    end
end
sig_factor = 0.2/mean(im_loc(im_loc(:)>0));
im_loc3d = reshape(im_loc,256,256,3);
figure; imagesc(im_loc3d*sig_factor); axis equal tight;



end
