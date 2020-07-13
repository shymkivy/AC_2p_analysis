function [data_flat, num_elem] = f_concat_dsets(dset, field_name, dim_cat)
% concatenates field data withing datasets  
% structure is dataset{1-n}.field -> concatenated field across datasets

f_names = fieldnames(dset{1});

if ~sum(strcmp(field_name, f_names))
    warning([field_name ' is not a field in the array']);
    return;
end

num_dsets = numel(dset);

if num_dsets == 1
    data_flat = dset{1}.(field_name);
    return;
end

num_indx = zeros(num_dsets,1);
siz_1 = size(dset{1}.(field_name));
num_indx(1) = siz_1(dim_cat);
siz_1(dim_cat) = [];

% check sizes of dsets
for n_dset = 2:num_dsets
    siz_n = size(dset{n_dset}.(field_name));
    num_indx(n_dset) = siz_n(dim_cat);
    siz_n(dim_cat) = [];
    if ~prod(siz_1 == siz_n)
        warning('other dimensions are non consistent in size');
        return;
    end
end

combined_dsets = cell(num_dsets,1);
for n_dset = 1:num_dsets
    combined_dsets{n_dset} = dset{n_dset}.(field_name);    
end

data_flat = cat(1,combined_dsets{:});
num_elem = sum(num_indx);

% if ~isempty(varargin)
%     f_names2 = varargin';
%     for n_arg = numel(f_names2):-1:1
%         if ~sum(strcmp(varargin{n_arg},f_names))
%             warning([f_names2{n_arg} ' input to f_combine_dsets() does not exist'])
%             f_names2(n_arg) = [];
%             
%         end
%     end
% else
%     f_names2 = f_names;
% end


% for n_field = 1:numel(f_names2)
%     
%     if strcmp(f_names2(n_field),'data_vec')
%         data_flat.(f_names2{n_field}) = zeros(siz_d, 'single');
%     elseif strcmp(f_names2(n_field),'trial_types')
%     end
%     
%     siz_d = size(dset.c_data{1}.(f_names2{n_field}));
%     inx_inx = (siz_d == dset.num_cells(1));
%     siz_d(inx_inx) = dset.num_cells_total;
%     data_flat.(f_names2{n_field}) = zeros(siz_d, 'single');
%     siz_d_start = ones(numel(siz_d),1);
%     siz_d_end = siz_d;
%     
%     if numel(siz_d) == 3
%         data_flat.(f_names2{n_field})(siz_d_start(1):siz_d_end(1),siz_d_start(2):siz_d_end(2),siz_d_start(3):siz_d_end(3)) = 1;
%     end
% 
% end


end