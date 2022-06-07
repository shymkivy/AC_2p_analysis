function Y = f_smooth_movie2(Y, smooth_std)
% fastest, uses GPU

siz = size(Y);
smooth_std1 = smooth_std;
siz1 = siz;

g = gpuDevice;
mem_use = g.FreeMemory/1e1;

Y = single(Y);
for n_sm = 1:3
    if smooth_std1(1)>0
        % make kernel
        sm_std1 = smooth_std1(1);
        kernel_half_size = ceil(sqrt(-log(0.05)*2*sm_std1^2));
        gaus_win = (-kernel_half_size:kernel_half_size)';
        gaus_kernel = exp(-((gaus_win).^2)/(2*sm_std1^2));
        gaus_kernel = gpuArray(gaus_kernel/sum(gaus_kernel));
        
        %figure; plot(gaus_kernel)
        
        norm_line = ones(siz1(1),1);
        norm_line_sm = gpuArray(conv(norm_line, gaus_kernel, 'same'));
        
        max_batch_size = floor(mem_use/siz1(1));
        num_batch = ceil((siz1(2)*siz1(3))/max_batch_size);
        batch_lims = ceil(linspace(0,siz1(2)*siz1(3),num_batch+1));

        Y = reshape(Y, siz1(1), siz1(2)*siz1(3));
        Y_batch = cell(num_batch,1);
        for n_b = 1:num_batch
            batch_data = gpuArray(Y(:,(batch_lims(n_b)+1):batch_lims(n_b+1)));
            batch_data = convn(batch_data, gaus_kernel, 'same')./norm_line_sm;
            Y_batch{n_b} = gather(batch_data);
        end
        Y = reshape(cat(2, Y_batch{:}), siz1(1), siz1(2), siz1(3));
    end
    smooth_std1 = smooth_std1([2 3 1]);
    siz1 = siz1([2 3 1]);
    Y = permute(Y, [2 3 1]);
end
clear batch_data;
Y = uint16(Y);

end