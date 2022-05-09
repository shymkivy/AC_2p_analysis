function image_out = f_mc_apply_frame_shift(image, shifts, params)

if ~exist('params', 'var')
    params = struct;
end

if ~isfield(params, 'zero_sides')
    params.zero_sides = 1;
end


[d1, d2, T] = size(image);

Lx = -d2/2:(d2/2-1);
Ly = -d1/2:(d1/2-1);

Lxi = ifftshift(Lx);
Lyi = ifftshift(Ly);

[Lxi2,Lyi2] = meshgrid(Lxi,Lyi);

image_out = zeros(d1,d2,T);

for n_t = 1:T
    image_ft = fft2(image(:,:,n_t));
    sh_corr = image_ft .* exp(-1i * 2*pi* (Lxi2*shifts(n_t,1)/d2 + Lyi2*shifts(n_t,2)/d1));
    image_out(:,:,n_t) = real(ifft2(sh_corr));
    if params.zero_sides
        if shifts(n_t,1) <= -1
            image_out(:,(end-floor(shifts(n_t,1))):end,n_t) = 0;
        elseif shifts(n_t,1) >= 1
            image_out(:,1:floor(shifts(n_t,1)),n_t) = 0;
        end
        if shifts(n_t,2) <= -1
            image_out((end-floor(shifts(n_t,2))):end,:,n_t) = 0;
        elseif shifts(n_t,2) >= 1
            image_out(1:floor(shifts(n_t,2)),:,n_t) = 0;
        end
    end
end

end