function f_mpl_plot_spatial_tunning(cdata, ops, cond_name, n_dset)
  
%[num_cells, n_dset] = max(cdata.num_cells);
num_cells = cdata.num_cells(n_dset);

for n_pl = 1:cdata.num_planes(n_dset)
    A = cdata.OA_data{n_dset,n_pl}.est.A(:,cdata.OA_data{n_dset,n_pl}.proc.comp_accepted);
    A = A/max(A(:));
    %A = reshape(A,256,256,num_cells);

    figure; imagesc(reshape(mean(A,2),256,256));
    axis equal tight;
    title(sprintf('%s %s %s,mpl%d, %d cells', cond_name,ops.paradigm_type, ops.file_names.(cond_name){n_dset},n_pl,num_cells), 'interpreter', 'none');


    %% freq tuning plots
    jet_map = jet(10);

    figure;

    rc_freq_onset = cdata.peak_tuned_trials_onset{n_dset}(:,1:10);
    rc_freq_onset_mag = cdata.tuning_all{n_dset}.peak_tuning_onset.fr_peak_mag_ave_z(:,1:10);

    subplot(2,3,1);
    if_plot_tuning_map(A,rc_freq_onset,rc_freq_onset_mag,jet_map);
    title(sprintf('Onset tuned cells n=%d', sum(logical(sum(rc_freq_onset,2)))), 'interpreter', 'none');

    rc_freq_offset = cdata.peak_tuned_trials_offset{n_dset}(:,1:10);
    rc_freq_offset_mag = cdata.tuning_all{n_dset}.peak_tuning_offset.fr_peak_mag_ave_z(:,1:10);

    subplot(2,3,2);
    if_plot_tuning_map(A,rc_freq_offset,rc_freq_offset_mag,jet_map);
    title(sprintf('Offset tuned cells n=%d', sum(logical(sum(rc_freq_offset,2)))), 'interpreter', 'none');

    rc_both = logical(rc_freq_onset + rc_freq_offset);
    rc_mag_both = (rc_freq_onset.*rc_freq_onset_mag + rc_freq_offset.*rc_freq_offset_mag);

    subplot(2,3,3);
    if_plot_tuning_map(A,rc_both,rc_mag_both,jet_map);
    title(sprintf('Onset and Offset tuned cells n=%d', sum(logical(sum(rc_both,2)))), 'interpreter', 'none');

    %figure; imagesc(1,ops.control_carrier_freq/1000, reshape(jet(100),100,1,3));

    %% deviance plot
    dd_map = [.8 .8 .8; .2 .6 1; 1 .2 .2];
    dd_map1 = [dd_map; dd_map];


    rc_mmn_onset = cdata.peak_tuned_trials_onset_ctx{n_dset};
    rc_mmn_onset_mag = cdata.tuning_all{n_dset}.peak_tuning_onset.fr_peak_mag_ave_z(:,cdata.ctx_mmn{n_dset});

    subplot(2,3,4);
    if_plot_tuning_map(A,rc_mmn_onset,rc_mmn_onset_mag,dd_map1);
    title(sprintf('Onset ctx tuned cells n=%d', sum(logical(sum(rc_mmn_onset,2)))), 'interpreter', 'none');

    rc_mmn_offset = cdata.peak_tuned_trials_offset_ctx{n_dset};
    rc_mmn_offset_mag = cdata.tuning_all{n_dset}.peak_tuning_offset.fr_peak_mag_ave_z(:,cdata.ctx_mmn{n_dset});

    subplot(2,3,5);
    if_plot_tuning_map(A,rc_mmn_offset,rc_mmn_offset_mag,dd_map1);
    title(sprintf('Offset ct tuned cells n=%d', sum(logical(sum(rc_mmn_offset,2)))), 'interpreter', 'none');

    rc_mmn_both = logical(rc_mmn_onset + rc_mmn_offset);
    rc_mmn_both_mag = rc_mmn_onset.*rc_mmn_onset_mag + rc_mmn_offset.*rc_mmn_offset_mag;
    subplot(2,3,6);
    if_plot_tuning_map(A,rc_mmn_both,rc_mmn_both_mag,dd_map1);
    title(sprintf('On and offset ctx tuned cells n=%d', sum(logical(sum(rc_mmn_both,2)))), 'interpreter', 'none');

    suptitle(sprintf('%s %s %s, mpl%d, n=%d cells total', cond_name,ops.paradigm_type, ops.file_names.(cond_name){n_dset},n_pl,num_cells));
end

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
imagesc(im_loc3d*sig_factor); axis equal tight;



end
