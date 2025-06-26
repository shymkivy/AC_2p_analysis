function f_dv_tuning_across_dsets(app)

if strcmpi(app.regdatatouseDropDown.Value, 'caiman reg')
    source = 'register_roi_caiman_load';
elseif strcmpi(app.regdatatouseDropDown.Value, 'gui reg')
    source = 'register_roi';
end

n_pl = app.mplSpinner.Value;

resp_type = app.ResposivecellstypeDropDown.Value;

data = app.data;
ddata = data(app.current_data_idx,:);
data2 = data(strcmpi(data.mouse_id, ddata.mouse_id),:);
data3 = data2(data2.FOV_num == ddata.FOV_num,:);

num_dsets = numel(data3.mouse_id);

all_combs = nchoosek(1:num_dsets,2);

for n_comb = 1:size(all_combs,1)
    dset1 = data3(all_combs(n_comb,1),:);
    dset2 = data3(all_combs(n_comb,2),:);

    reg_data1 = dset1.(source){n_pl};
    reg_data2 = dset2.(source){n_pl};
    
    if ~isempty(reg_data1)
    
        reg_mat = [reg_data1.reg_cell_idx, reg_data2.reg_cell_idx];
        keep_cells = ~(sum(isnan(reg_mat),2));
        reg_mat2 = reg_mat(keep_cells,:);
        num_cells = size(reg_mat2,1);
        
        stats1 = dset1.stats{n_pl};
        stats2 = dset2.stats{n_pl};

        if strcmpi(resp_type, 'Peaks')
            resp_vals1 = stats1.peak_vals;
            resp_vals2 = stats2.peak_vals;
        elseif strcmpi(resp_type, 'Onset')
            resp_vals1 = stats1.onset_vals;
            resp_vals2 = stats2.onset_vals;
        elseif strcmpi(resp_type, 'Offset')
            resp_vals1 = stats1.offset_vals;
            resp_vals2 = stats2.offset_vals;
        end
        tuning1 = resp_vals1(:,1:10);
        tuning2 = resp_vals2(:,1:10);

        tuning1_cut = tuning1(reg_mat2(:,1),:);
        tuning2_cut = tuning2(reg_mat2(:,2),:);
        
        figure;
        plot(max(tuning1_cut, [], 2), max(tuning2_cut, [], 2), '.');
        title(sprintf('tuning max mag comparison; %s', resp_type));
        xlabel('dset 1');
        ylabel('dset 2');

        figure;
        subplot(1,2,1);
        imagesc(tuning1_cut);
        subplot(1,2,2);
        imagesc(tuning2_cut);
        sgtitle(sprintf('freq tuning; %', resp_type));

        similarity = 1 - diag(pdist2(tuning1_cut, tuning2_cut, 'correlation'));
        
        tuning_mag = mean([max(tuning1_cut, [], 2), max(tuning2_cut, [], 2)], 2);

        figure;
        plot(tuning_mag, similarity, '.');
        xlabel('mean max tuning mag');
        ylabel('correlation');
        title(sprintf('correlation vs tuning mag; %s', resp_type));
        
        est1 = dset1.OA_data{n_pl}.est;
        est2 = dset2.OA_data{n_pl}.est;

        A1 = reshape(full(est1.A(:,reg_mat2(:,1))), [est1.dims(1), est1.dims(2), num_cells]);
        A2 = reshape(full(est2.A(:,reg_mat2(:,2))), [est2.dims(1), est2.dims(2), num_cells]);
        
        A1_im = sum(A1,3)';
        A2_im = sum(A2,3)';
        norm_fac = max(max([A1_im, A2_im]))/1.5;

        rgb_im = cat(3, A1_im/norm_fac, A2_im/norm_fac, A1_im/norm_fac);

        figure; imagesc(rgb_im); 
        title(sprintf('%s im%d_%s(m) vs im%d_%s(g), pl%d; %d cells; %s', ddata.mouse_id{1}, dset1.im_num, dset1.dset_name{1}, dset2.im_num , dset2.dset_name{1}, n_pl, num_cells, app.regdatatouseDropDown.Value), 'interpreter', 'none')
    else
        fprintf('no data\n')
    end
end

end