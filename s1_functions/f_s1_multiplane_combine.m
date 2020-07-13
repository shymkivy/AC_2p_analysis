function trace = f_s1_multiplane_combine(mpl_struct)

num_planes = numel(mpl_struct);
num_frames_total = 0;
for n_pl = 1:num_planes
    sizS = size(mpl_struct{n_pl});
    num_frames_total = num_frames_total + max(sizS);
end
num_chan = min(sizS);
if num_chan == sizS(1)
    trace = zeros(num_chan, num_frames_total);
elseif num_chan == sizS(2)
    trace = zeros(num_frames_total,num_chan);
end
for n_pl = 1:num_planes
    if num_chan == sizS(1)
        trace(:, n_pl:num_planes:num_frames_total) = mpl_struct{n_pl};
    elseif num_chan == sizS(2)
        trace(n_pl:num_planes:num_frames_total,:) = mpl_struct{n_pl};
    end
end
end


