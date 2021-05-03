function [shifted_scaled_trace] = align_volt_by_scale_shift2(trace, scaling_factor, shift)
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
    scaled_trace = interp1(x,trace,xq);
    if num_chan == 1
        scaled_trace = scaled_trace';
    end
elseif scaling_factor == 1
    scaled_trace = trace;
end

% now adjust for the shift
if shift < 0
    % remove voltage data
    shifted_scaled_trace = scaled_trace(1-shift:end,:);
elseif shift > 0
    % insert voltage data
    pad_val = scaled_trace(1);
    temp_padding = ones(shift,num_chan)*pad_val;
    if d1<d2
        shifted_scaled_trace = [temp_padding', scaled_trace];
    else
        shifted_scaled_trace = [temp_padding; scaled_trace];
    end
elseif shift == 0
    shifted_scaled_trace = scaled_trace;
end


end