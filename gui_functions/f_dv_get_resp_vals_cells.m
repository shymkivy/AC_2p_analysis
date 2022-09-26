function [resp_vals, resp_cells] = f_dv_get_resp_vals_cells(stats, tn_all, feature_type, resp_cell_select, resp_thr)


if ~exist('resp_cell_select', 'var')
    resp_cell_select = 'Resp marg'; % all, resp marg, resp split
end

if ~exist('feature_type', 'var')
    feature_type = 'Peaks';
end

num_trials = numel(tn_all);
num_cells = sum(stats.num_cells);

num_slice = 1;
if strcmpi(feature_type, 'Peaks')
    vals1 = stats.peak_val_all.*stats.peak_in_resp_win;
    resp_thr2 = mean(stats.resp_thresh_peak(:,stats.lim_win_frames),2)*resp_thr;
elseif strcmpi(feature_type, 'Onset')
    vals1 = stats.onset_vals;
    resp_thr2 = stats.resp_thresh_onset*resp_thr;
elseif strcmpi(feature_type, 'Offset')
    vals1 = stats.offset_vals;
    resp_thr2 = stats.resp_thresh_offset*resp_thr;
elseif strcmpi(feature_type, 'OnOff')
    num_slice = 2;
    % takes whichever is more significant by  z score
    vals1 = cat(3,stats.onset_vals, stats.offset_vals);
    resp_thr2 = [stats.resp_thresh_onset, stats.resp_thresh_offset];
else
    error('feature type undefined');
end

resp_cells1 = cell(num_trials, 1);
vals = cell(num_trials, 1);

for n_tn = 1:num_trials
    tn1 = tn_all(n_tn);
    
    idx_out = logical(sum(vals1(:,tn1,1) > resp_thr2(:,1),2));
    vals_out = vals1(:,tn1,1);
    
    if num_slice > 1
        vals_z1 = vals1(:,tn1,1)./resp_thr2(:,1);
        vals_z2 = vals1(:,tn1,2)./resp_thr2(:,2);
        
        vals2 = vals1(:,tn1,2);
        idx2 = logical(sum(vals1(:,tn1,2) > resp_thr2(:,2),2));
        
        [~, idx3] = max([vals_z1, vals_z2], [], 2);
        
        idx_out(idx3==2) = idx2(idx3==2);
        vals_out(idx3==2) = vals2(idx3==2);
    end

    resp_cells1{n_tn} = idx_out;
    vals{n_tn} = vals_out;
end

if strcmpi(resp_cell_select, 'All')
    resp_cells = true(num_cells, num_trials);
elseif strcmpi(resp_cell_select, 'Resp marg')
    resp_cells = repmat(logical(sum(cat(2,resp_cells1{:}),2)), [1, num_trials]);
elseif strcmpi(resp_cell_select, 'Resp split')
    resp_cells = cat(2,resp_cells1{:});
end

resp_vals = cell(num_trials, 1);
for n_tn = 1:num_trials
    resp_vals{n_tn} = vals{n_tn}(resp_cells(:,n_tn));
end


end