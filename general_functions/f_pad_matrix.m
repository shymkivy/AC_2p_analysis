function mat_out = f_pad_matrix(mat_in, num_pad, rep)

if ~exist('num_pad', 'var')
    num_pad = 1;
end

if ~exist('rep', 'var')
    rep = 0;
end

siz = size(mat_in);

if rep
    mat_out1 = [ones(siz(1), num_pad, siz(3)).*mat_in(:,1,:), mat_in, ones(siz(1), num_pad, siz(3)).*mat_in(:,end,:)];
    mat_out = [ones(num_pad, siz(2)+2*num_pad, siz(3)).*mat_out1(1,:,:); mat_out1; ones(num_pad, siz(2)+2*num_pad, siz(3)).*mat_out1(end,:,:)];
else
    mat_out1 = [zeros(siz(1), num_pad, siz(3)), mat_in, zeros(siz(1), num_pad, siz(3))];
    mat_out = [zeros(num_pad, siz(2)+2*num_pad, siz(3)); mat_out1; zeros(num_pad, siz(2)+2*num_pad, siz(3))];
end

end