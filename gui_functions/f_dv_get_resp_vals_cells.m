function [resp_vals, resp_cells] = f_dv_get_resp_vals_cells(app, stats, tn_all, feature_type, resp_cell_select, resp_thr)

if ~exist('feature_type', 'var') || isempty(feature_type)
    feature_type = app.ResposivecellstypeDropDown.Value;
end

if ~exist('resp_cell_select', 'var') || isempty(resp_cell_select)
    resp_cell_select = app.ResponsivecellsselectDropDown.Value; % all, resp marg, resp split
end

if ~exist('resp_thr', 'var') || isempty(resp_thr)
    resp_thr = app.RespthreshEditField.Value; % all, resp marg, resp split
end

num_trials = numel(tn_all);
num_cells = sum([stats.num_cells]);

num_slice = 1;

if strcmpi(feature_type, 'Peaks')
    vals1 = cat(1, stats.peak_vals).*cat(1, stats.peak_in_resp_win);
    resp_thr2 = cat(1, stats.peak_resp_thresh)*resp_thr;
elseif strcmpi(feature_type, 'Onset')
    vals1 = cat(1, stats.onset_vals);
    resp_thr2 = cat(1, stats.onset_resp_thresh)*resp_thr;
elseif strcmpi(feature_type, 'Offset')
    vals1 = cat(1, stats.offset_vals);
    resp_thr2 = cat(1, stats.offset_resp_thresh)*resp_thr;
elseif strcmpi(feature_type, 'OnOff')
    num_slice = 2;
    % takes whichever is more significant by  z score
    vals1 = cat(3,cat(1, stats.onset_vals), cat(1, stats.offset_vals));
    resp_thr2 = [cat(1, stats.onset_resp_thresh), cat(1, stats.offset_resp_thresh)];
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