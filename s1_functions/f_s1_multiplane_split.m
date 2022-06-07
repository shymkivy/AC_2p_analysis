function [mpl_stack, mpl_num_frames] = f_s1_multiplane_split(trace, num_planes)

SizT = size(trace);
num_frames = max(SizT);
num_frames = floor(num_frames/num_planes)*num_planes;

mpl_stack = cell(1,num_planes);
mpl_num_frames = zeros(1,num_planes);

for n_pl = 1:num_planes
    mpl_num_frames(n_pl) = ceil((num_frames-n_pl+1)/num_planes);
    if SizT(1) > SizT(2)
        mpl_stack{n_pl} = trace(n_pl:num_planes:num_frames,:);
    elseif SizT(1) < SizT(2)
        mpl_stack{n_pl} = trace(:,n_pl:num_planes:num_frames);
    end
end

end