function shift_xy = f_mc_compute_frame_shift(target, image, params)
% F(target) = F(image) * exp(-i * 2pi * shift);

if ~exist('params', 'var')
    params = struct;
end

if ~isfield(params, 'shift_method')
    params.shift_method = 'phase_angle';
end

[d1, d2, T] = size(image);

im_2d = reshape(image, [d1*d2, T]);
im_mean = mean(im_2d);
im_std = std(im_2d);
image_n = (image - im_mean)/im_std;

targ_mean = mean(target(:));
targ_std = std(target(:));
target_n = (target - targ_mean)/targ_std;

target_ft = ifftshift(fft2(target_n));

shift_xy = zeros(T, 2);
if strcmpi(params.shift_method, 'phase_angle')
    for n_t = 1:T
        image_ft = fft2(image_n(:,:,n_t));
        corr1 = target_ft./image_ft;
        [Fx, Fy] = gradient(real(log(corr1)/(-1i)/2/pi));
        shift_xy(n_t, :) = [median(Fx(:))*d2, median(Fy(:))*d1];
    end
else
    for n_t = 1:T
        image_ft = ifftshift(fft2(image_n(:,:,n_t)));
        im_conv = ifft2(exp(-1i*imag(target_ft))*exp(-1i*imag(image_ft)));
        [Fx, Fy] = gradient(real(log(corr1)/(-1i)/2/pi));
        shift_xy(n_t, :) = [median(Fx(:))*d2, median(Fy(:))*d1];
    end
end

end