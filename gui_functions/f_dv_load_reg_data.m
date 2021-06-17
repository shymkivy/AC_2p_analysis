function f_dv_load_reg_data(app)

disp('Loading reg data...');
data1 = load(app.regdatapathEditField.Value);
reg_data = data1.data_all;

num_dset = numel(reg_data);
num_regions = size(reg_data(1).wf_mapping_regions_coords,2);

for n_dset = 1:num_dset
    region_means = zeros(num_regions, 2); 
    coords = reg_data(n_dset).wf_mapping_regions_coords;
    for n_reg = 1:num_regions
        temp1 = coords(:,n_reg);
        region_means(n_reg, :) = mean(cat(1,temp1{:}));
    end
    reg_data(n_dset).wf_region_means = region_means;
end

app.reg_data = reg_data;
disp('Done');

end