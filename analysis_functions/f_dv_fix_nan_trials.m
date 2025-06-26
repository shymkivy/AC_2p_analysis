function [data_out, rem_idx] = f_dv_fix_nan_trials(data_in, method)

data_out = data_in;
num_cells = size(data_out,1);
rem_idx = false(num_cells,1);

hasnan1 = logical(sum(isnan(data_in),2));
if strcmpi(method, 'randn')
    for n_cell = 1:num_cells
        if hasnan1(n_cell)
            data3 = data_in(n_cell,:);
            nan_idx = isnan(data3);
            vals1 = randn(sum(nan_idx),1);
            vals1 = vals1 * std(data3(~nan_idx));
            vals1 = vals1 + mean(data3(~nan_idx));
            data_out(n_cell,nan_idx) = vals1;
        end
    end
elseif strcmpi(method, 'mean')
    for n_cell = 1:num_cells
        if hasnan1(n_cell)
            data3 = data_in(n_cell,:);
            nan_idx = isnan(data3);
            data_out(n_cell,nan_idx) = mean(data3(~nan_idx));
        end
    end
elseif strcmpi(method, 'remove')
    data_out = data_in(~hasnan1,:);
    rem_idx = hasnan1;
end

end