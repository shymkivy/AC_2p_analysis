function f_dv_low_dim_distances(app)

dist_metric = app.LDdistmethodDropDown.Value; % pca isomap

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

num_dsets = size(data,1);

num_cont = 8;

tn_all = [1:num_cont 18 19 20 28 29 30]; %f_dv_get_trial_number(app);

title_tag1 = sprintf('%s; dist %s', title_tag, dist_metric);


data_all = cell(num_dsets,1);
cell_dset_idx = cell(num_dsets,1);
num_cells_all = zeros(num_dsets,1);
for n_dset = 1:num_dsets
    stats1 = cat(1,data(n_dset,:).stats{n_pl});

    if strcmpi(app.ResponsivecellsselectDropDown.Value, 'All')
        resp_cell_sel = 'All';
    else
        resp_cell_sel = 'Resp marg';
    end

    [~, resp_vals] = f_dv_get_resp_vals_cells(app, stats1, tn_all, [], resp_cell_sel);

    resp_vals2 = cat(2,resp_vals{:});
    num_cells_all(n_dset) = size(resp_vals2,1);
    cell_dset_idx{n_dset} = ones(num_cells_all(n_dset),1)*n_dset;
   
    data_all{n_dset} = resp_vals2;
end


data_all2 = cat(1,data_all{:});
hasnan1 = logical(sum(isnan(data_all2),2));
data_all2 = data_all2(~hasnan1,:);

cell_dset_idx2 = cat(1,cell_dset_idx{:});
MMN_freq = data.MMN_freq;

tn_d = [18, 19, 20; 28, 29, 30];
MMMN_id = [2 2 2; 1 1 1];
MMMN_id_opp = [1 1 1; 2 2 2];

[num_mmn, num_tnd] = size(tn_d);
dist_all = cell(num_dsets, num_mmn, num_tnd);

for n_dset = 1:num_dsets
    title_tag2 = sprintf('%s; dset%d; %dcells', title_tag1, data.idx(n_dset), num_cells_all(n_dset));
    resp_vals = data_all{n_dset};

    cont_vec = resp_vals(:,1:num_cont);
    
    for n_mmn = 1:num_mmn
        for n_tnd = 1:num_tnd
            tn_idx = tn_all == tn_d(n_mmn, n_tnd);
            vec1 = resp_vals(:,tn_idx);
            dist_all{n_dset, n_mmn, n_tnd} = pdist2(cont_vec', vec1', dist_metric);
        end
    end
end

dist_all_same = cell(num_dsets, n_mmn, n_tnd);
dist_all_opp = cell(num_dsets, n_mmn, n_tnd);

pad = 2;
for n_dset = 1:num_dsets
    
    MMN_freq1 = MMN_freq{n_dset};
    for n_mmn = 1:num_mmn
        for n_tnd = 1:num_tnd
            MMN_freq2 = MMN_freq1(MMMN_id(n_mmn, n_tnd));
            MMN_freq2_opp = MMN_freq1(MMMN_id_opp(n_mmn, n_tnd));
            dist_temp = dist_all{n_dset, n_mmn, n_tnd};

            dist_temp2 = ones(pad*2+1, 1);
            for n_pad = -pad:pad
                if and(MMN_freq2+n_pad >= 1 , MMN_freq2+n_pad <= num_cont)
                    dist_temp2(pad+1+n_pad) = dist_temp(MMN_freq2+n_pad);
                end
            end
            dist_all_same{n_dset, n_mmn, n_tnd} = dist_temp2;

            
            dist_temp2 = ones(pad*2+1, 1);
            for n_pad = -pad:pad
                if and(MMN_freq2_opp+n_pad >= 1 , MMN_freq2_opp+n_pad <= num_cont)
                    dist_temp2(pad+1+n_pad) = dist_temp(MMN_freq2_opp+n_pad);
                end
            end
            dist_all_opp{n_dset, n_mmn, n_tnd} = dist_temp2;
        end
    end
end



x_lab = -pad:pad;
figure; hold on;
for n_tnd = 1:num_tnd
    dist_temp = cat(2, dist_all_same{:, :, n_tnd});
    dist_temp_opp = cat(2, dist_all_opp{:, :, n_tnd});
    
    plot(x_lab, dist_temp, '-', 'color', [app.ops.context_types_all_colors2{tn_d(1, n_tnd)} 0.2])
    plot(x_lab, mean(dist_temp,2), 'o-', 'color', app.ops.context_types_all_colors2{tn_d(1, n_tnd)}, 'linewidth', 2)

    plot(x_lab, dist_temp_opp, '--', 'color', [app.ops.context_types_all_colors2{tn_d(1, n_tnd)}, 0.2])
    plot(x_lab, mean(dist_temp_opp,2), 'o--', 'color', app.ops.context_types_all_colors2{tn_d(1, n_tnd)}, 'linewidth', 2)
end
xlabel('Frequency around MMN');
ylabel([dist_metric ' distance']);
title(title_tag1)


%%
MMN_freq5 = cat(1, MMN_freq{:});

% first location in mmn is 2, second is 1
high_idx = 2 - double(MMN_freq5(:,1) < MMN_freq5(:,2));
low_idx = 3 - high_idx;

dist_temp_high_same = cell(num_dsets, n_tnd);
dist_temp_high_opp = cell(num_dsets, n_tnd);
dist_temp_low_same = cell(num_dsets, n_tnd);
dist_temp_low_opp = cell(num_dsets, n_tnd);

for n_dset = 1:num_dsets
    for n_tnd = 1:num_tnd
        dist_temp_high_same{n_dset, n_tnd} = dist_all_same{n_dset, high_idx(n_dset), n_tnd};
        dist_temp_high_opp{n_dset, n_tnd} = dist_all_opp{n_dset, high_idx(n_dset), n_tnd};

        dist_temp_low_same{n_dset, n_tnd} = dist_all_same{n_dset, low_idx(n_dset), n_tnd};
        dist_temp_low_opp{n_dset, n_tnd} = dist_all_opp{n_dset, low_idx(n_dset), n_tnd};
    end
end


x_lab = -pad:pad;
figure; hold on;
for n_tnd = 1:num_tnd

    dist_temp = cat(2, dist_temp_high_same{:, n_tnd});
    dist_temp_opp = cat(2, dist_temp_high_opp{:, n_tnd});
    
    plot(x_lab, dist_temp, '-', 'color', [app.ops.context_types_all_colors2{tn_d(1, n_tnd)} 0.2])
    plot(x_lab, mean(dist_temp,2), 'o-', 'color', app.ops.context_types_all_colors2{tn_d(1, n_tnd)}, 'linewidth', 2)

    plot(x_lab, dist_temp_opp, '--', 'color', [app.ops.context_types_all_colors2{tn_d(1, n_tnd)}, 0.2])
    plot(x_lab, mean(dist_temp_opp,2), 'o--', 'color', app.ops.context_types_all_colors2{tn_d(1, n_tnd)}, 'linewidth', 2)
end
xlabel('High frequency MMN');
ylabel([dist_metric ' distance']);
title([title_tag1 '; high MMN'])

x_lab = -pad:pad;
figure; hold on;
for n_tnd = 1:num_tnd

    dist_temp = cat(2, dist_temp_low_same{:, n_tnd});
    dist_temp_opp = cat(2, dist_temp_low_opp{:, n_tnd});
    
    plot(x_lab, dist_temp, '-', 'color', [app.ops.context_types_all_colors2{tn_d(1, n_tnd)} 0.2])
    plot(x_lab, mean(dist_temp,2), 'o-', 'color', app.ops.context_types_all_colors2{tn_d(1, n_tnd)}, 'linewidth', 2)

    plot(x_lab, dist_temp_opp, '--', 'color', [app.ops.context_types_all_colors2{tn_d(1, n_tnd)}, 0.2])
    plot(x_lab, mean(dist_temp_opp,2), 'o--', 'color', app.ops.context_types_all_colors2{tn_d(1, n_tnd)}, 'linewidth', 2)
end
xlabel('Low frequency MMN');
ylabel([dist_metric ' distance']);
title([title_tag1 '; low MMN'])

% 
% for n_dset = 1:num_dsets
%     figure; hold on;
%     for n_tnd = 1:num_tnd
%         if tn_d(n_tnd) > 20
%             pl_line = 'o--';
%         else
%             pl_line = 'o-';
%         end
%         plot(dist_all{n_dset, n_tnd}, pl_line, 'color', app.ops.context_types_all_colors2{tn_d(n_tnd)})
%     end
%     title(title_tag2, 'interpreter', 'none');
% end

end