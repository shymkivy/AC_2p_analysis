function str_out = f_array_to_str(array_in)

num_el = numel(array_in);
cell1 = cell(2, num_el);

for n_el = 1:num_el
    cell1{1, n_el} = num2str(array_in(n_el));
    cell1{2, n_el} = ', ';
end

cell2 = cell1(:)';
cell2(end) = [];

str_out = cat(2,cell2{:});

end