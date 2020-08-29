function smooth_data = f_smooth_gauss2(data, sigma_frames, normalize)

if ~exist('sigma_frames', 'var') || isempty(sigma_frames)
    sigma_frames = 1; % default
end

if ~exist('normalize', 'var') || isempty(normalize)
    normalize = 0; % default
end

% make kernel
kernel_half_size = ceil(sqrt(-log(0.05)*2*sigma_frames^2));
gaus_win = -kernel_half_size:kernel_half_size;
gaus_kernel = exp(-((gaus_win).^2)/(2*sigma_frames^2));
gaus_kernel = gaus_kernel/sum(gaus_kernel);

%figure; plot(gaus_win, gaus_kernel);

% extract the number of cells and the number of frames
num_frames = size(data,2);
num_cells = size(data,1);

smooth_data=zeros(num_cells, num_frames);
for n_cell=1:num_cells    % size(data,1)
    temp_data = data(n_cell,:);

    % convolve gaussian
    temp_data = conv2(temp_data, gaus_kernel, 'same');

    % normalize
    if normalize
        temp_data = temp_data/max(temp_data);
    end
    
    smooth_data(n_cell,:) = temp_data;
end

end