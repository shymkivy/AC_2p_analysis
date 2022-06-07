function shift_xy = f_mc_compute_frame_shift(target, image, params)
% F(target) = F(image) * exp(-i * 2pi * shift);

if ~exist('params', 'var')
    params = struct;
end

if ~isfield(params, 'shift_method')
    params.shift_method = 'phase_angle'; % phase_angle
end

[d1, d2, T] = size(image);
d1f = d1*2-1;
d1_cent = ceil(d1f/2);

im_2d = reshape(image, [d1*d2, T]);
im_mean = mean(im_2d);
im_std = std(im_2d);
image_n = (image - im_mean)/im_std;

targ_mean = mean(target(:));
targ_std = std(target(:));
target_n = (target - targ_mean)/targ_std;

target_ft = fft2(rot90(target_n,2), d1f, d1f);

shift_xy = zeros(T, 2);
if strcmpi(params.shift_method, 'phase_angle')
    for n_t = 1:T
        image_ft = fft2(rot90(image_n(:,:,n_t),2), d1f, d1f);
        corr1 = target_ft./image_ft;
        [Fx, Fy] = gradient(angle(ifftshift(corr1)));
        %[Fx, Fy] = gradient(real(log(corr1./abs(corr1))/(-1i)/2/pi));
        shift_xy(n_t, :) = [median(Fx(:))*d1f, median(Fy(:))*d1f];
    end
else
    for n_t = 1:T
        image_ft = fft2(image_n(:,:,n_t), d1f, d1f);
        im_conv1 = ifft2(target_ft.*image_ft);
        [~, idx1] = max(im_conv1(:));
        [row, col] = ind2sub([d1f, d1f], idx1);
        shift_xy(n_t, :)= d1_cent - [col, row];
    end
end

end