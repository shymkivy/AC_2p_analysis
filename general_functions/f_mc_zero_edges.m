function Y = f_mc_zero_edges(Y, dsall1_use, ds_base_all)
% y-x orientations
num_frames = size(dsall1_use,1);

dsall1_use_r = round(dsall1_use + ds_base_all);

for n_fr = 1:num_frames          
    if dsall1_use_r(n_fr,1) < 0
        Y(1:-dsall1_use_r(n_fr,1),:,n_fr) = 0;
    elseif dsall1_use_r(n_fr,1) > 0
        Y((end-dsall1_use_r(n_fr,1)):end,:,n_fr) = 0;
    end
    if dsall1_use_r(n_fr,2) < 0
        Y(:,1:-dsall1_use_r(n_fr,2),n_fr) = 0;
    elseif dsall1_use_r(n_fr,2) > 0
        Y(:,(end-dsall1_use_r(n_fr,2)):end,n_fr) = 0;
    end
end

end