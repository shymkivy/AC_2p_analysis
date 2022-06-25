function AC_data2 = f_s0_parse_tab_data(params)

AC_data = readtable(params.dset_table_fpath);

%%

AC_data = AC_data(~isnan(AC_data.idx),:);

AC_data = AC_data(AC_data.do_proc == 1,:);

%%
limit1 = params.limit;

fields1 = fields(limit1);

AC_data2 = AC_data;
for n_fl = 1:numel(fields1)
    field_cur = fields1{n_fl};
    if isfield(limit1, field_cur)
        tag1 = limit1.(field_cur);
        if numel(tag1)
            AC_data2 = AC_data2(strcmpi(AC_data2.(field_cur), tag1),:);
        end
    end
end


end