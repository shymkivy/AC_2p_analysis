function [resp_cells, resp_vals, resp_vals_full, resp_locs] = f_dv_get_resp_vals_cells(app, stats, tn_all, feature_type, resp_cell_select, resp_thr)

if ~exist('feature_type', 'var') || isempty(feature_type)
    feature_type = app.ResposivecellstypeDropDown.Value;
end

if ~exist('resp_cell_select', 'var') || isempty(resp_cell_select)
    resp_cell_select = app.ResponsivecellsselectDropDown.Value; % all, resp marg, resp split
end

if ~exist('resp_thr', 'var') || isempty(resp_thr)
    resp_thr = app.RespthreshEditField.Value; % all, resp marg, resp split
end

num_cells = sum([stats.num_cells]);
num_trials = numel(tn_all);

num_slice = 1;

if strcmpi(feature_type, 'Peaks')
    % get vals and replace what is outside of window with means
    vals1 = cat(1, stats.peak_vals);
    in_win_idx = cat(1, stats.peak_in_resp_win);
    tr_ave_mean = squeeze(cat(1, stats.trial_ave_lim_win_mean));
    vals1(~in_win_idx) = tr_ave_mean(~in_win_idx);
    resp_thr2 = cat(1, stats.peak_resp_thresh)*resp_thr;
    locs = cat(1, stats.peak_loc);
    locs(~in_win_idx) = NaN;
elseif strcmpi(feature_type, 'Onset')
    vals1 = cat(1, stats.onset_vals);
    resp_thr2 = cat(1, stats.onset_resp_thresh)*resp_thr;
    locs = stats.onset_loc;
elseif strcmpi(feature_type, 'Offset')
    vals1 = cat(1, stats.offset_vals);
    resp_thr2 = cat(1, stats.offset_resp_thresh)*resp_thr;
    locs = stats.offset_loc;
elseif strcmpi(feature_type, 'OnOff')
    num_slice = 2;
    % takes whichever is more significant by  z score
    vals1 = cat(3,cat(1, stats.onset_vals), cat(1, stats.offset_vals));
    resp_thr2 = [cat(1, stats.onset_resp_thresh), cat(1, stats.offset_resp_thresh)];
    locs = ones(size(vals1));
    locs(:,:,1) = stats.onset_loc;
    locs(:,:,2) = stats.offset_loc;
else
    error('feature type undefined');
end

resp_cells1 = cell(num_trials, 1);
vals = cell(num_trials, 1);
locs2 = cell(num_trials, 1);
for n_tn = 1:num_trials
    tn1 = tn_all(n_tn);
    
    idx_out = logical(sum(vals1(:,tn1,1) > resp_thr2(:,1),2));
    vals_out = vals1(:,tn1,1);
    locs_out = locs(:,tn1,1);
    
    if num_slice > 1
        vals_z1 = vals1(:,tn1,1)./resp_thr2(:,1);
        vals_z2 = vals1(:,tn1,2)./resp_thr2(:,2);
        
        vals2 = vals1(:,tn1,2);
        idx2 = logical(sum(vals1(:,tn1,2) > resp_thr2(:,2),2));
        locs2 = locs(:,:,2);
        
        [~, idx3] = max([vals_z1, vals_z2], [], 2);
        
        idx_out(idx3==2) = idx2(idx3==2);
        vals_out(idx3==2) = vals2(idx3==2);
        locs_out(idx3==2) = locs2(idx3==2);
    end

    resp_cells1{n_tn} = idx_out;
    vals{n_tn} = vals_out;
    locs2{n_tn} = locs_out;
end

if strcmpi(resp_cell_select, 'All')
    resp_cells = true(num_cells, num_trials);
elseif strcmpi(resp_cell_select, 'Resp marg')
    resp_cells = repmat(logical(sum(cat(2,resp_cells1{:}),2)), [1, num_trials]);
elseif strcmpi(resp_cell_select, 'Resp split')
    resp_cells = cat(2,resp_cells1{:});
end

resp_vals_full = cat(2,vals{:});

resp_locs = cat(2,locs2{:});

resp_vals = cell(num_trials, 1);
for n_tn = 1:num_trials
    resp_vals{n_tn} = vals{n_tn}(resp_cells(:,n_tn));
end

end