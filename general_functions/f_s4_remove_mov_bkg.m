function Y_alln = f_s4_remove_mov_bkg(Y_all, method)

[d1, d2, T] = size(Y_all);

Y_all = single(reshape(Y_all, d1*d2, T));

if method == 1
    % either simple mean removal from each axis
    Y_alln = Y_all - mean(Y_all,1);
    Y_alln = Y_alln - mean(Y_alln,2);
elseif method == 2
    % or compute a rank 1 background component with some iterating
    num_it = 4;

    Ya2d_spat = mean(Y_all,2);
    Ya2d_spat = Ya2d_spat/norm(Ya2d_spat);
    Ya2d_temp = mean(Y_all,1);
    Ya2d_temp = Ya2d_temp/norm(Ya2d_temp);

    Ya2d_spat_curr = Ya2d_spat;
    Ya2d_temp_curr = Ya2d_temp;
    Ya2d_spat2 = cell(num_it,1);
    Ya2d_temp2 = cell(num_it,1);
    err1_spat = zeros(num_it,1);
    err1_temp = zeros(num_it,1);
    for n_it = 1:num_it
        Ya2d_temp_curr = (Ya2d_spat_curr' * Y_all);
        Ya2d_temp_curr = Ya2d_temp_curr/norm(Ya2d_temp_curr);

        Ya2d_spat_curr = (Y_all * Ya2d_temp_curr');
        Ya2d_spat_curr = Ya2d_spat_curr/norm(Ya2d_spat_curr);

        Ya2d_temp2{n_it} = Ya2d_temp_curr;
        Ya2d_spat2{n_it} = Ya2d_spat_curr;

        err1_spat(n_it) = mean(Ya2d_spat - Ya2d_spat2{n_it});
        err1_temp(n_it) = mean(Ya2d_temp - Ya2d_temp2{n_it});
    end

    figure;
    subplot(2,1,1);
    plot(err1_spat); title('spatial mean diff')
    subplot(2,1,2);
    plot(err1_temp); title('temporal mean diff')

    Ya2d_temp3 = Ya2d_temp2{end};
    Ya2d_spat3 = (Y_all * Ya2d_temp3');

    Y_bkg = Ya2d_spat3*Ya2d_temp3;
    Y_alln = Y_all - Y_bkg;
    clear Y_bkg;
end

Y_all = reshape(Y_all, d1, d2, T);
Y_alln = reshape(Y_alln, d1, d2, T);

end