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


%% compute all transformations

interp_k = 0;%2; % number of times splits the coords into two (k)
% if used ad increase fac, 185*2^interp_r_fac - 2^interp_r_fac +1
% num points in between 2^k - 1

%%
all_mice = unique(app.data.mouse_tag, 'stable');

%%
%anchor_dset = app.anchordsetSpinner.Value;
anchor_mouse_tag = '10_2_18'; %all_mice{anchor_dset};
%%
anch_idx = strcmpi({reg_data.mouse_tag}, anchor_mouse_tag);
anchor_reg = reg_data(:,anch_idx);
%%
region_means_all = cat(3,reg_data.wf_region_means);
anchor_means = region_means_all(:,:,anch_idx);

%%
mouse_tforms = cell(numel(all_mice),1);
for n_ms = 1:numel(all_mice)
    mouse_tag2 = all_mice{n_ms};

    %%
    current_rdata_idx = strcmpi(mouse_tag2,{reg_data.mouse_tag});
    current_means = region_means_all(:,:,current_rdata_idx);
%     current_means_in = region_means_all_in(:,:,current_rdata_idx);
    
    current_tform_wf = fitgeotrans(current_means,anchor_means ,'nonreflectivesimilarity');
    current_tform_wf.T(3,1:2) = current_tform_wf.T(3,1:2)*2^interp_k - 2^interp_k + 1;
    
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
        
        cdata = f_dv_compute_cdata(app, params);
        
        %% load images and tform
        tform = rdata2.regions_tforms.tform;

        %% load the contours and get coords
        %accepted_cells = mdata2.stats{1}{n_pl}.accepted_cells;
        A = mdata.OA_data{1}.est.A(:,cdata.accepted_cells);

        %num_cells = mdata2.stats{1}.num_cells; % n_pl
        num_cells = cdata.num_cells;
        coords = ones(num_cells,3);
        [~, ind1] = max(A);
        [coords(:,1), coords(:,2)] = ind2sub([256 256], ind1);

        %% interpolate

        tform_in = tform;
        tform_in.T(3,1:2) = tform_in.T(3,1:2)*2^interp_k - 2^interp_k + 1;

        tform_in.T = tform_in.T * current_tform_wf.T;

        coords_in = coords;
        coords_in(:,1:2) = coords_in(:,1:2)*2^interp_k - 2^interp_k + 1;

        %% register 
        coords_tf = (coords_in*tform_in.T);
        
        app.data(n_dset,:).registered_data{1}.coords = coords_tf;
    end
end
disp('Done');

%%
num_freq = 3;
im_wf = anchor_reg.wf_mapping_im{num_freq};


cell_labels = cell(num_dsets,1);
for n_dset = 1:num_dsets
    if ~isempty(app.data(n_dset,:).registered_data{1})
        coords_tf = app.data(n_dset,:).registered_data{1}.coords;
        num_cells = size(coords_tf,1);
        pos_min_dist = zeros(num_cells,4);
        for n_cell = 1:num_cells
            for n_reg = 1:4
                pos1 = anchor_reg.region_borders{n_reg}.Position;
                pos_min_dist(n_cell, n_reg) = min(sqrt(sum((pos1 - coords_tf(n_cell,1:2)).^2,2)));
            end
        end
        [~, reg_idx] = min(pos_min_dist,[],2);
        cell_labels{n_dset} = reg_idx;
        app.data(n_dset,:).registered_data{1}.reg_labels = reg_idx;
    end
end

figure; imagesc(im_wf); axis equal tight; hold on;
for n_reg = 1:4
    pos1 = anchor_reg.region_borders{n_reg}.Position;
    plot(pos1(:,1), pos1(:,2),'Color', app.ops.cond_colors{n_reg}, 'LineWidth', 2);
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

end