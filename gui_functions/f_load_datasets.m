function [data, ops] = f_load_datasets(ops)

disp('Loading data...');
data = ops.AC_data;

% gather the file names of specified data

if isfield(ops, 'num_dsets_load')
    num_dsets = ops.num_dsets_load;
else
    num_dsets = numel(data.area);
end
if ops.waitbar
    if isfield(ops, 'app')
        app = ops.app;
    else
        app = [];
    end
    wb = f_waitbar_initialize(app, 'Loading data...');
end

has_data = false(num_dsets,1);
for n_dset = 1:num_dsets
    ddata = data(n_dset,:);   
    fname_tag = ddata.dset_name{1};
    date_tag = ddata.mouse_tag{1};
    
    % load proc data
    temp_proc = dir([ops.file_dir, ['\' '*' fname_tag '*' date_tag '*' ops.processed_data_tag '.mat']]);
    if ~isempty(temp_proc)
        temp_load = load([ops.file_dir '\' temp_proc.name]);
        data.proc_data{n_dset} = temp_load.data;
        data.proc_ops{n_dset} = temp_load.ops;
    end
    

    % load OA
    temp_OA = dir([ops.file_dir, ['\' '*' fname_tag '*' date_tag '*' ops.OA_output_tag '.mat']]);
    if ~isempty(temp_OA)
        has_data(n_dset) = 1;
        data.num_planes(n_dset) = numel(temp_OA);
        for n_pl = 1:data.num_planes(n_dset)
            if data.num_planes(n_dset) > 1
                if isempty(strfind(temp_OA(n_pl).name, ['pl'  num2str(n_pl)]))
                    error('File name and planes dont match');
                end
            end
            temp_data1 = load([ops.file_dir '\' temp_OA(n_pl).name]);
            
            % for throw bad components for memory effic
            idx_acc = logical(temp_data1.proc.comp_accepted);
            idx_acc(temp_data1.proc.idx_manual_bad) = 0;
            idx_acc(temp_data1.proc.idx_manual) = 1;
            temp_data1.est = if_cut_bad_comp(temp_data1.est, idx_acc);
            temp_data1.proc = if_cut_bad_comp(temp_data1.proc, idx_acc);
            temp_data1.proc.idx_cell_acc = idx_acc;
            temp_data1.proc.deconv.smooth_dfdt = if_cut_bad_comp(temp_data1.proc.deconv.smooth_dfdt, idx_acc);
            if isfield(temp_data1.proc.deconv, 'df_f')
                temp_data1.proc.deconv.df_f = if_cut_bad_comp(temp_data1.proc.deconv.df_f, idx_acc);
            end
            if isfield(temp_data1.proc.deconv, 'c_foopsi')
                temp_data1.proc.deconv.c_foopsi = if_cut_bad_comp(temp_data1.proc.deconv.c_foopsi, idx_acc);
            end
            if isfield(temp_data1.proc.deconv, 'MCMC')
                temp_data1.proc.deconv.MCMC = if_cut_bad_comp(temp_data1.proc.deconv.MCMC, idx_acc);
            end
                
            
            data.OA_data{n_dset, n_pl} = temp_data1;
        end
    end
    
    if ops.waitbar
        f_waitbar_update(wb, n_dset/num_dsets, sprintf('Loading %d/%d',n_dset,num_dsets));
    end
end

if ops.waitbar
    f_waitbar_close(wb);
end
data = data(has_data,:);

% more stuff

max_planes = max(data.num_planes);
num_dsets = size(data,1);

for n_var = 1:numel(ops.save_var_list_pl)
    var1 = ops.save_var_list_pl{n_var};
    if ~sum(strcmpi(data.Properties.VariableNames, var1))
        data.(var1) = cell(num_dsets,max_planes);
    end
end

for n_var = 1:numel(ops.save_var_list)
    var1 = ops.save_var_list{n_var};
    if ~sum(strcmpi(data.Properties.VariableNames, var1))
        data.(var1) = cell(num_dsets,1);
    end
end

data.registered_data = cell(num_dsets,max_planes);

end

%%
function struct_out = if_cut_bad_comp(struct_in, comp_idx)

struct_out = struct_in;

num_comp = numel(comp_idx);
            
fields1 = fields(struct_in);
for n_fl = 1:numel(fields1)
    field_temp = fields1{n_fl};
    siz1 = size(struct_in.(field_temp));
    if siz1(1) == num_comp
        struct_out.(field_temp) = struct_in.(field_temp)(comp_idx,:);
    elseif siz1(2) == num_comp
        struct_out.(field_temp) = struct_in.(field_temp)(:,comp_idx);
    end
end


end