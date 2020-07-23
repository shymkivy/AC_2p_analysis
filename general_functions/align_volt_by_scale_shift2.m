function [shifted_scaled_dat_proc] = align_volt_by_scale_shift2(trace, scaling_factor, shift)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This fuction takes in a voltage trace and shift and scale factos 
%
%   Last update: 11/27/19 Yuriy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% align here

shift = round(shift);

[d1, d2] = size(trace);
num_samp = max([d1, d2]);
num_chan = min([d1, d2]);

if scaling_factor > 1 || scaling_factor < 1
    % insert missing
    x = ((1:num_samp)-1)*scaling_factor+1;
    xq = 1:1:num_samp*scaling_factor;
    scaled_dat_proc = interp1(x,trace,xq);

elseif scaling_factor == 1
    scaled_dat_proc = trace;
end


% now adjust for the shift
if shift < 0
    % remove voltage data
    shifted_scaled_dat_proc = scaled_dat_proc(1-shift:end,:);
elseif shift > 0
    % insert voltage data
    pad_val = scaled_dat_proc(1);
    temp_padding = ones(shift,num_chan)*pad_val;
    shifted_scaled_dat_proc = [temp_padding; scaled_dat_proc];
elseif shift == 0
    shifted_scaled_dat_proc = scaled_dat_proc;
end


end