function Y = f_smooth_movie_double(Y, smooth_std)
% fast but needs ram

siz = size(Y);
smooth_std1 = smooth_std;
siz1 = siz;

for n_sm = 1:3
    if smooth_std1(1)>0
        % make kernel
        sm_std1 = smooth_std1(1);
        kernel_half_size = ceil(sqrt(-log(0.05)*2*sm_std1^2));
        gaus_win = (-kernel_half_size:kernel_half_size)';
        gaus_kernel = exp(-((gaus_win).^2)/(2*sm_std1^2));
        gaus_kernel = gaus_kernel/sum(gaus_kernel);
        
        %figure; plot(gaus_kernel)
        
        norm_line = ones(siz1(1),1);
        norm_line_sm = conv(norm_line, gaus_kernel, 'same');
 
        Y = reshape(Y, siz1(1), siz1(2)*siz1(3));
        Y = convn(Y, gaus_kernel, 'same')./norm_line_sm;
        Y = reshape(Y, siz1(1), siz1(2), siz1(3));
        
    end
    smooth_std1 = smooth_std1([2 3 1]);
    siz1 = siz1([2 3 1]);
    Y = permute(Y, [2 3 1]);
end

end