function Y = f_mc_cut_high_vals(Y, thresh, per_frame)

if ~exist('thresh', 'var')
    thresh = 0;
end

if ~exist('per_frame', 'var')
    per_frame = 1;
end

[d1, d2, T] = size(Y);
Y = reshape(Y, [d1*d2, T]);

if thresh>0
    if per_frame
        valcut = prctile(Y, (1-thresh)*100, 1);
        for n_t = 1:T
            idx1 = Y(:,n_t)>valcut(n_t);
            Y(idx1, n_t) = valcut(n_t);
        end
    else
        valcut = prctile(Y(:), (1-thresh)*100);
        for n_t = 1:T
            idx1 = Y(:,n_t)>valcut;
            Y(idx1, n_t) = valcut;
        end
    end
end

Y = reshape(Y, [d1, d2, T]);

end