function [shifted_scaled_dat_proc] = align_volt_by_scale_shift(dat_proc, scaling_factor, shift)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This fuction takes in a voltage trace with coordinates from two traces
%   that need to be scaled and shifted. [beg1, end1; beg2, end2]. This
%   scales the trace with coordinates 2 to coordinates 1.
%
%   Last update: 3/16/18 Yuriy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% align here

shift = round(shift);

temp_scaling_index_array = 1:size(dat_proc,1);

scaled_index_array = round(temp_scaling_index_array*scaling_factor);


if scaling_factor > 1
    % insert missing

    % find the locations of missing frames
    missing_frames = find(diff(scaled_index_array) == 2);

    % create new larger array
    scaled_dat_proc = zeros(length(scaled_index_array)+length(missing_frames),size(dat_proc,2));

    % create array with locations that need t obe duplicated
    scaled_dat_proc_index = zeros(length(scaled_index_array)+length(missing_frames),1);
    scaled_dat_proc_index(missing_frames+(1:length(missing_frames))) = 1;

    % fill scaled array with data
    scaled_dat_proc(logical(scaled_dat_proc_index),:) = dat_proc(missing_frames,:);
    scaled_dat_proc(logical(1-scaled_dat_proc_index),:) = dat_proc;


elseif scaling_factor < 1

    % find the locations with duplicated frames
    extra_frames = diff(scaled_index_array) == 0;

    scaled_dat_proc = dat_proc;

    scaled_dat_proc(extra_frames,:) = [];
elseif scaling_factor == 1
    
    scaled_dat_proc = dat_proc;
end


% now adjust for the shift

if shift < 0
    % remove voltage data
    shifted_scaled_dat_proc = scaled_dat_proc(1-shift:end,:);
elseif shift > 0
    % insert voltage data
    shifted_scaled_dat_proc = [zeros(shift,size(scaled_dat_proc,2)); scaled_dat_proc];
elseif shift == 0
    shifted_scaled_dat_proc = scaled_dat_proc;
end


end