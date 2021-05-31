function [smooth_dfdt_data, sig_std] = f_smooth_dfdt(data, params)
% this calculates derivative of dfof of calcium traces, with Jordans method
% yuriy 11/15/19

%% check inputs
if ~exist('params', 'var')
    params = struct;
end

if isfield(params, 'do_smooth')
    do_smooth = params.do_smooth;
else
    do_smooth = 1; % default
end

if isfield(params, 'sigma_frames')
    sigma_frames = params.sigma_frames;
else
    sigma_frames = 1; % default
end

if isfield(params, 'rectify')
    rectify = params.rectify;
else
    rectify = 1; % default
end

if isfield(params, 'normalize')
    normalize = params.normalize;
else
    normalize = 0; % default
end    

if isfield(params, 'threshold')
    threshold = params.threshold;
else
    threshold = 0; % default in Z
end

%% process
% make kernel
kernel_half_size = ceil(sqrt(-log(0.05)*2*sigma_frames^2));
gaus_win = -kernel_half_size:kernel_half_size;
gaus_kernel = exp(-((gaus_win).^2)/(2*sigma_frames^2));
gaus_kernel = gaus_kernel/sum(gaus_kernel);

% extract the number of cells and the number of frames
num_frames = size(data,2);
num_cells = size(data,1);

smooth_dfdt_data=zeros(num_cells, num_frames);
sig_std = zeros(num_cells,1);
for n_cell=1:num_cells    % size(data,1)
    temp_data = data(n_cell,:);
    
    % derivative
    temp_data = [0 diff(temp_data)];

    % convolve gaussian
    if do_smooth
        temp_data = conv2(temp_data, gaus_kernel, 'same');
    end
    
    % normalize
    if normalize
        temp_data = temp_data/max(temp_data);
    end
    
    % rectify
    if rectify
        temp_data = max(temp_data,0);
    end
    
    % threshold
    if threshold
        % first get the z factor
        temp_thresh_data = temp_data;
        temp_thresh_data(temp_thresh_data <= 0) = [];
        sig_std(n_cell) = rms(temp_thresh_data);
        temp_data(temp_data<=threshold*sig_std(n_cell)) = threshold*sig_std(n_cell);
        temp_data = temp_data - threshold*sig_std(n_cell);
    end
    
    smooth_dfdt_data(n_cell,:) = temp_data;
end


end