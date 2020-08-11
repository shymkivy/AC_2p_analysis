function param_out = f_get_param(struct, param_name)

if isfield(struct, param_name)
    param_out = struct.(param_name);
else
    param_out = [];
end

end