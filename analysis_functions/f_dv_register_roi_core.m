function reg_out = f_dv_register_roi_core(data, throw_thresh)

%throw_thresh = 0.5;

num_dsets = numel(data.mouse_id);
num_planes = data.num_planes(1);

%%
A_flat_all = cell(num_planes, num_dsets);
A_all3 = cell(num_planes, num_dsets);
A_all2 = cell(num_planes, num_dsets);
A_all2_sp = cell(num_planes, num_dsets);
num_cells = zeros(num_planes, num_dsets);
for n_dset = 1:num_dsets
    for n_pl = 1:num_planes
        data3 = data(n_dset,:);
        est = data3.OA_data{n_pl}.est;
        proc = data3.OA_data{n_pl}.proc;
        num_cells(n_pl, n_dset) = sum(proc.comp_accepted);
        A1 = est.A(:,logical(proc.comp_accepted));
        A_full = full(A1);
        A_flat = sum(A_full,3);
        %figure; imagesc(A_flat);
        %title(num2str(fov_idx(n_idx)));
        A_all2_sp{n_pl, n_dset} = sparse(A1);
        A_all2{n_pl, n_dset} = A_full;
        A_all3{n_pl, n_dset} = reshape(A_full, est.dims(1), est.dims(2), []);
        A_flat_all{n_pl, n_dset} = A_flat/max(A_flat(:))*2;
    end
end

%% compute matches for all dset pairs
all_pairs = nchoosek(1:num_dsets, 2);
num_pairs = size(all_pairs,1);
match_idx_all = cell(num_planes,num_pairs);
match_scores_all = cell(num_planes,num_pairs);
for n_pl = 1:num_planes
    for n_pair = 1:size(all_pairs,1)
        dset1 = all_pairs(n_pair, 1);
        dset2 = all_pairs(n_pair, 2);
        
        A1_full = A_all2_sp{n_pl, dset1};
        A2_full = A_all2_sp{n_pl, dset2};
        
        num_cells1 = num_cells(n_pl, dset1);
        num_cells2 = num_cells(n_pl, dset2);
        match_idx = [(1:num_cells1)', nan(num_cells1,1)];
        match_scores = zeros(num_cells1,1);
        for n_cell1 = 1:num_cells1
            match_scores2 = sum(A1_full(:,n_cell1).*A2_full);
            [max_score, n_cell2] = max(match_scores2);
            if max_score > throw_thresh
                match_idx(n_cell1,2) =  n_cell2;
                match_scores(n_cell1) = full(max_score);
            end
        end
        num_rem = 0;
        not_reg_idx = false(num_cells2,1);
        for n_cell2 = 1:num_cells2
            idx2 = n_cell2 == match_idx(:,2);
            idx2_locs = find(idx2);
            % remove duplicate matches
            if sum(idx2) > 1
                [~, idx3] = max(match_scores(idx2_locs));
                rem_locs = idx2_locs;
                rem_locs(idx3) = [];
                match_idx(rem_locs,2) = nan;
                match_scores(rem_locs) = 0;
                num_rem = num_rem + numel(rem_locs);
            end
            if ~sum(idx2)
                not_reg_idx(n_cell2) = 1;
            end
        end
        rem_mat = [nan(sum(not_reg_idx),1), find(not_reg_idx)];
        match_idx_all{n_pl, n_pair} = [match_idx; rem_mat];
        match_scores_all{n_pl, n_pair} = [match_scores; zeros(sum(not_reg_idx),1)];
    end
end

%% combine pairs across all dsets with greedy method
match_idx_mpl = cell(num_planes,1);
match_scores_mpl = cell(num_planes,1);

for n_pl = 1:num_planes

    match_idx_stack = cell(num_pairs,1);
    match_scores_stack = cell(num_pairs,1);
    all_pairs_stack = cell(num_pairs,1);
    for n_pair = 1:size(all_pairs,1)
        
        idx1 = logical(match_scores_all{n_pl, n_pair});
        match_idx_stack{n_pair} = match_idx_all{n_pl, n_pair}(idx1,:);
        match_scores_stack{n_pair} = match_scores_all{n_pl, n_pair}(idx1);
        all_pairs_stack{n_pair} = repmat(all_pairs(n_pair,:), [sum(idx1), 1]);
    end
    
    match_idx_stack2 = cat(1, match_idx_stack{:});
    match_scores_stack2 = cat(1, match_scores_stack{:});
    all_pairs_stack2 = cat(1, all_pairs_stack{:});
    
    [~, idx1] = sort(match_scores_stack2, 'descend');
    
    match_idx_stack3 = match_idx_stack2(idx1,:);
    match_scores_stack3 = match_scores_stack2(idx1,:);
    all_pairs_stack3 = all_pairs_stack2(idx1,:);
    
    num_match = numel(idx1);

    n_match = 1;
    
    new_row = nan(1,num_dsets);
    new_score = zeros(1,num_dsets);
    new_row(all_pairs_stack3(n_match, 1)) = match_idx_stack3(n_match, 1);
    new_row(all_pairs_stack3(n_match, 2)) = match_idx_stack3(n_match, 2);
    new_score(all_pairs_stack3(n_match, 2)) = match_scores_stack3(n_match);
    match_idx_allall = new_row;
    match_score_allall = new_score;
    
    for n_match = 2:num_match
        added = 0;
        idx1 = match_idx_allall(:,all_pairs_stack3(n_match, 1)) == match_idx_stack3(n_match, 1);
        idx2 = match_idx_allall(:,all_pairs_stack3(n_match, 2)) == match_idx_stack3(n_match, 2);
        if sum(idx1) && ~sum(idx2)
            % if first cell already exists and not second
            if isnan(match_idx_allall(idx1,all_pairs_stack3(n_match, 2)))
                % if second space is nan, fill
                match_idx_allall(idx1,all_pairs_stack3(n_match, 2)) = match_idx_stack3(n_match, 2);
                match_score_allall(idx1,all_pairs_stack3(n_match, 2)) = max(match_scores_stack3(n_match), match_score_allall(idx1,all_pairs_stack3(n_match, 2)));
                added = 1;
            else
                % there is conflict, skipping
                %fprintf('skipping missmatch dset%d to dset%d, cell%d to cell%d, pl%d\n', all_pairs_stack3(n_match, 1), all_pairs_stack3(n_match, 2), match_idx_stack3(n_match, 1), match_idx_stack3(n_match, 2), n_pl);
                added = 1;
            end
        elseif sum(idx2) && ~sum(idx1)
            % if second cell already exists and not first
            if isnan(match_idx_allall(idx2,all_pairs_stack3(n_match, 1)))
                % if first is nan, fill
                match_idx_allall(idx2,all_pairs_stack3(n_match, 1)) = match_idx_stack3(n_match, 1);
                match_score_allall(idx2,all_pairs_stack3(n_match, 2)) = max(match_scores_stack3(n_match), match_score_allall(idx2,all_pairs_stack3(n_match, 2)));
                added = 1;
            else
                % there is conflict, skipping
                %fprintf('skipping missmatch dset%d to dset%d, cell%d to cell%d, pl%d\n', all_pairs_stack3(n_match, 1), all_pairs_stack3(n_match, 2), match_idx_stack3(n_match, 1), match_idx_stack3(n_match, 2), n_pl);
                added = 1;
            end
        elseif sum(idx1) && sum(idx2)
            % if both cells already exist
            if ~sum(idx1.*idx2)
                % if first and second are in different locations
                % look for conflicts
                temp_match_idx = match_idx_allall(logical(idx1 + idx2),:);
                conf_all = false(num_dsets,1);
                for n_dset = 1:num_dsets
                    if ~sum(isnan(temp_match_idx(:,n_dset)))
                        if temp_match_idx(1,n_dset) ~= temp_match_idx(2,n_dset)
                            conf_all(n_dset) = 1;
                        end
                    end
                end
                % if no conflicts, merge rows, else leave as is
                if ~sum(conf_all)
                    for n_dset = 1:num_dsets
                        if ~isnan(match_idx_allall(idx2, n_dset))
                            % copy vals
                            match_idx_allall(idx1, n_dset) = match_idx_allall(idx2, n_dset);
                            match_score_allall(idx1, n_dset) = match_score_allall(idx2, n_dset);
                        end
                    end
                    match_score_allall(idx1,all_pairs_stack3(n_match, 2)) = max(match_scores_stack3(n_match), match_score_allall(idx1,all_pairs_stack3(n_match, 2)));
                    % remove the second
                    match_idx_allall(idx2,:) = [];
                    match_score_allall(idx2,:) = [];
                    added = 1;
                else
                    %fprintf('skipping missmatch dset%d to dset%d, cell%d to cell%d, pl%d\n', all_pairs_stack3(n_match, 1), all_pairs_stack3(n_match, 2), match_idx_stack3(n_match, 1), match_idx_stack3(n_match, 2), n_pl);
                    added = 1;
                end
            else
                % already filled
                match_score_allall(idx1,all_pairs_stack3(n_match, 2)) = max(match_scores_stack3(n_match), match_score_allall(idx1,all_pairs_stack3(n_match, 2)));
                added = 1;
            end
        end
        
        % add new rows for new matches
        if ~added
            new_row = nan(1,num_dsets);
            new_score = zeros(1,num_dsets);
            new_row(all_pairs_stack3(n_match, 1)) = match_idx_stack3(n_match, 1);
            new_row(all_pairs_stack3(n_match, 2)) = match_idx_stack3(n_match, 2);
            new_score(all_pairs_stack3(n_match, 2)) = match_scores_stack3(n_match);
            match_idx_allall = [match_idx_allall; new_row];
            match_score_allall = [match_score_allall; new_score];
        end
    end
    
    % fill missing cells with no reg
    for n_dset = 1:num_dsets
        dset_cells_all = (1:num_cells(n_pl, n_dset))';
        dset_cells = match_idx_allall(:,n_dset);
        dset_cells2 = dset_cells(~isnan(dset_cells));
        if numel(unique(dset_cells2)) ~= numel(dset_cells2)
            %fprintf('duplicates in dset %d pl%d', n_dset, n_pl);
        end
        
        idx1 = sum(dset_cells_all == dset_cells2',2);
        missing_cells = dset_cells_all(~idx1);
        match_idx_extra = nan(numel(missing_cells), num_dsets);
        match_idx_extra(:,n_dset) = missing_cells;
        match_idx_allall = [match_idx_allall; match_idx_extra];
    end
    
    % sort 
    for n_dset = 1:num_dsets
        [~, idx1] = sort(match_idx_allall(:,(num_dsets - n_dset + 1)));
        match_idx_allall = match_idx_allall(idx1,:);
    end
    
    match_idx_mpl{n_pl} = match_idx_allall;
	match_scores_mpl{n_pl} = match_score_allall;
end



%% combine paris across all dsets combine then fix method (not good yet, has dupl)
% 
% match_idx_mpl2 = cell(num_planes,1);
% match_scores_mpl2 = cell(num_planes,1);
% for n_pl = 1:num_planes
%     match_idx_allall = [match_idx_all{n_pl, 1}, nan(size(match_idx_all{n_pl, 1},1), num_dsets-2)];
%     match_score_allall = [zeros(size(match_idx_all{n_pl, 1},1), 1), match_scores_all{n_pl, 1}, zeros(size(match_idx_all{n_pl, 1},1), num_dsets-2)];
%     
%     for n_pair = 2:size(all_pairs,1)
%         dset1 = all_pairs(n_pair,1);
%         dset2 = all_pairs(n_pair,2);
%         
%         temp_data1 = match_idx_all{n_pl, n_pair};
%         temp_score1 = match_scores_all{n_pl, n_pair};
%         for n_cell = 1:max(temp_data1(:,1))
%             
%             idx_all = match_idx_allall(:,dset1) == n_cell;
%             idx1 = temp_data1(:,1) == n_cell;
%             if sum(idx_all)
%                 if ~isnan(temp_data1(idx1,2))            
%                     if isnan(match_idx_allall(idx_all,dset2))
%                         match_idx_allall(idx_all,dset2) = temp_data1(idx1,2);
%                         match_score_allall(idx_all,dset2) = max(temp_score1(idx1), match_score_allall(idx_all,dset2));
%                         %fprintf('new val dset%d vs dset%d, cell%d\n', dset1, dset2, n_cell)
%                     else
%                         if ~(match_idx_allall(idx_all,dset2) == temp_data1(idx1,2))
%                             fprintf('missmatch dset%d vs dset%d, cell%d, pl%d\n', dset1, dset2, n_cell, n_pl);
%                             % here resolve missmatched alignments
%                             
% %                             temp_scores2 = match_score_allall(idx_all,:);
% %                             temp_idx2 = match_idx_allall(idx_all,:);
% %                             
% %                             old_row = nan(1, num_dsets);
% %                             old_row(1) = temp_idx2(1);
% %                             new_row = nan(1, num_dsets);
% %                             
% %                             [max_val, max_idx] = max(temp_scores2);
% %                             
% %                             if max_val > temp_score1(idx1)
% %                                 old_row(max_idx) = temp_idx2(max_idx);
% %                                 temp_scores2(max_idx) = 0;
% %                             else
% %                                 new_row(max_idx) = temp_idx2(max_idx);
% %                                 temp_scores2(max_idx) = 0;
% %                             end
%                             
%                             if (temp_score1(idx1) > match_score_allall(idx_all,dset1))
%                                 % if new connection stronger than old point
%                                 
%                                 new_row = nan(1, num_dsets);
%                                 new_row(dset1) = match_idx_allall(idx_all,dset1);
%                                 new_row(dset2) = temp_data1(idx1,2);
%                                 new_scores = zeros(1, num_dsets);
%                                 new_scores(dset2) = temp_score1(idx1);
%                                 match_idx_allall(idx_all,dset1) = nan;
%                                 match_score_allall(idx_all,dset1) = 0;
%                                 match_idx_allall = [match_idx_allall; new_row];
%                                 match_score_allall = [match_score_allall; new_scores];
%                             end
% 
%                         else
%                             %fprintf('same val dset%d vs dset%d, cell%d\n', dset1, dset2, n_cell)
%                         end
%                     end
%                 end
%             else
%                 %fprintf('val missing dset%d vs dset%d, cell%d\n', dset1, dset2, n_cell)
%                 new_row = nan(1, num_dsets);
%                 new_row(dset1) = temp_data1(idx1,1);
%                 new_row(dset2) = temp_data1(idx1,2);
%                 match_idx_allall = [match_idx_allall; new_row];
%                 new_scores = zeros(1, num_dsets);
%                 match_score_allall = [match_score_allall; new_scores];
%             end
%         end
%   
%     end
%     
%     n_dset = num_dsets;
%     for n_cell = 1:num_cells(n_pl, n_dset)
%         
%         idx_all = match_idx_allall(:,n_dset) == n_cell;
%         if ~sum(idx_all)
%             %fprintf('adding new row dset%d, cell%d\n', n_dset, n_cell)
%             new_row = nan(1, num_dsets);
%             new_row(n_dset) = n_cell;
%             match_idx_allall = [match_idx_allall; new_row];
%             new_scores = zeros(1, num_dsets);
%             match_score_allall = [match_score_allall; new_scores];
%         end
%         
%     end
%     
%     match_idx_mpl2{n_pl} = match_idx_allall;
%     match_scores_mpl2{n_pl} = match_score_allall;
% end

%%
reg_out = cell(num_dsets, num_planes);
for n_pl = 1:num_planes
    reg_data = struct();
    reg_data.fov_cell_idx = (1:size(match_idx_mpl{n_pl},1))';
    for n_dset = 1:num_dsets
        reg_data2 = reg_data;
        reg_data2.reg_cell_idx = match_idx_mpl{n_pl}(:,n_dset);
        reg_out{n_dset, n_pl} = reg_data2;
    end
end

end