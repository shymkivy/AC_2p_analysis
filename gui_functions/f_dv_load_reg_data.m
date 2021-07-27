function f_dv_load_reg_data(app)

disp('Loading reg data...');
data1 = load(app.regdatapathEditField.Value);
reg_data = data1.data_all;

num_reg = numel(reg_data);
num_regions = size(reg_data(1).wf_mapping_regions_coords,2);

for n_dset_reg = 1:num_reg
    region_means = zeros(num_regions, 2); 
    coords = reg_data(n_dset_reg).wf_mapping_regions_coords;
    for n_reg = 1:num_regions
        temp1 = coords(:,n_reg);
        region_means(n_reg, :) = mean(cat(1,temp1{:}));
    end
    reg_data(n_dset_reg).wf_region_means = region_means;
end

app.reg_data = reg_data;

%%
all_mice = unique(app.data.mouse_tag, 'stable');
%%
anchor_dset = app.anchordsetSpinner.Value;
anchor_mouse_tag = all_mice{anchor_dset};
%%
anch_idx = strcmpi({reg_data.mouse_tag}, anchor_mouse_tag);
anchor_reg = reg_data(:,anch_idx);
%%
region_means_all = cat(3,reg_data.wf_region_means);
anchor_means = region_means_all(:,:,anch_idx);

%%
borders_dset = 1;
borders_reg = app.reg_data(:,borders_dset);
borders = borders_reg.region_borders;
borders_means = region_means_all(:,:,borders_dset);

borders_tform_wf = fitgeotrans(borders_means,anchor_means ,'nonreflectivesimilarity');

borders_tf = borders;
for n_reg = 1:4
    bor_pos = borders_tf{n_reg}.Position(:,1:2);
    bor_pos2 = [bor_pos, ones(size(bor_pos,1),1)]*borders_tform_wf.T; 
    borders_tf{n_reg}.Position = bor_pos2(:,1:2);
end

%%
mouse_tforms = cell(numel(all_mice),1);
for n_ms = 1:numel(all_mice)
    mouse_tag2 = all_mice{n_ms};

    %%
    current_rdata_idx = strcmpi(mouse_tag2,{reg_data.mouse_tag});
    current_means = region_means_all(:,:,current_rdata_idx);
%     current_means_in = region_means_all_in(:,:,current_rdata_idx);
    
    current_tform_wf = fitgeotrans(current_means,anchor_means ,'nonreflectivesimilarity');
    
    mouse_tforms{n_ms} = current_tform_wf;
end

params = f_dv_gather_params(app);

num_dsets = size(app.data,1);
disp('Computing cell locations')
for n_dset = 1:num_dsets
    params.n_dset = n_dset;
    mdata = app.data(n_dset,:);
    
    current_rdata_idx = strcmpi(mdata.mouse_tag,{reg_data.mouse_tag});
    current_rdata = reg_data(current_rdata_idx);
    rdata = current_rdata.regions;
    
    current_tform_wf = mouse_tforms{current_rdata_idx};

    idx_r = strcmpi(mdata.area, [rdata.region_name]);
    rdata2 = rdata(:,idx_r);

    if ~isempty(rdata2.regions_tforms)
        %% load images and tform
        tform = rdata2.regions_tforms.tform;
        tform.T = tform.T * current_tform_wf.T;
        
        %% load the contours and get coords
        %accepted_cells = mdata2.stats{1}{n_pl}.accepted_cells;
        A = mdata.OA_data{1}.est.A(:,mdata.stats{1}.accepted_cells);

        %num_cells = mdata2.stats{1}.num_cells; % n_pl
        num_cells = mdata.stats{1}.num_cells;
        coords = ones(num_cells,3);
        [~, ind1] = max(A);
        [coords(:,1), coords(:,2)] = ind2sub([256 256], ind1);

        %% register 
        coords_tf = (coords*tform.T);
        
        app.data(n_dset,:).registered_data{1}.coords = coords_tf;
    end
end
disp('Done');

%%

% upsample borders to fill gaps
max_dist = 1;
borders_fix = borders;
for n_reg = 1:4
    pos1 = borders{n_reg}.Position;
    n_pos = 1;
    while n_pos < size(pos1,1)
        dist1 = sqrt(sum((pos1(n_pos+1,:) - pos1(n_pos,:)).^2));
        if dist1 > max_dist
            x_lin = linspace(pos1(n_pos,1), pos1(n_pos+1,1), ceil(dist1/max_dist)+1);
            y_lin = linspace(pos1(n_pos,2), pos1(n_pos+1,2), ceil(dist1/max_dist)+1);
            pos1 = [pos1(1:n_pos,:); [x_lin', y_lin']; pos1((n_pos+1):end,:)];
        end
        n_pos = n_pos + 1;
    end
    dist1 = sqrt(sum((pos1(end,:) - pos1(1,:)).^2));
    x_lin = linspace(pos1(end,1), pos1(1,1), ceil(dist1/max_dist)+1);
    y_lin = linspace(pos1(end,2), pos1(1,2), ceil(dist1/max_dist)+1);
    borders_fix{n_reg}.Position = [pos1; [x_lin', y_lin']];
end
app.border_coords = borders_fix;

% sort cells to bordered regions
cell_labels = cell(num_dsets,1);
for n_dset = 1:num_dsets
    if ~isempty(app.data(n_dset,:).registered_data{1})
        coords_tf = app.data(n_dset,:).registered_data{1}.coords;
        num_cells = size(coords_tf,1);
        pos_min_dist = zeros(num_cells,4);
        for n_cell = 1:num_cells
            for n_reg = 1:4
                pos1 = borders_fix{n_reg}.Position;
                pos_min_dist(n_cell, n_reg) = min(sqrt(sum((pos1 - coords_tf(n_cell,1:2)).^2,2)));
            end
        end
        [~, reg_idx] = min(pos_min_dist,[],2);
        cell_labels{n_dset} = reg_idx;
        app.data(n_dset,:).registered_data{1}.reg_labels = reg_idx;
    end
end

%num_freq = 3;
%im_wf = anchor_reg.wf_mapping_im{num_freq};
f1 = figure; %imagesc(im_wf); 
axis equal tight; hold on;
for n_reg = 1:4
    pos1 = borders_fix{n_reg}.Position;  
    plot(pos1(:,1), pos1(:,2), '.', 'Color', app.ops.cond_colors{n_reg}, 'LineWidth', 2);
end
for n_dset = 1:num_dsets
    if ~isempty(app.data(n_dset,:).registered_data{1})
        coords_tf = app.data(n_dset,:).registered_data{1}.coords;
        for n_reg = 1:4
            reg_idx = cell_labels{n_dset} == n_reg;
            plot(coords_tf(reg_idx,1),coords_tf(reg_idx,2), '.', 'color',app.ops.cond_colors{n_reg});
        end
    end
end
f1.Children.YDir = 'reverse';

f_dv_load_dset_from_data(app);
end