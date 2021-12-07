function Y_out = f_bidi_apply_shift(Y_in, shifts)

[d1, d2, T] = size(Y_in);
idx_1 = 1:2:d1;
idx_2 = 2:2:d1;

Y_out = Y_in;

for n_fr = 1:T
    frame_in = Y_in(:,:,n_fr);
    frame_out = frame_in - frame_in;

    rem1 = rem(shifts(n_fr),2);
    half_shift = (shifts(n_fr)-rem1)/2;

    % shift odd lines right
    shift1 = -half_shift;
    start_in1 = max([1-shift1 1]);
    end_in1 = min([d2-shift1, d2]);
    start_out1 = max([1+shift1 1]);
    end_out1 = min([d2+shift1, d2]);
    frame_out(idx_1,start_out1:end_out1) = frame_in(idx_1,start_in1:end_in1);


    % shift even lines left
    shift2 = half_shift + rem1;
    start_in2 = max([1-shift2 1]);
    end_in2 = min([d2-shift2, d2]);
    start_out2 = max([1+shift2 1]);
    end_out2 = min([d2+shift2, d2]);
    frame_out(idx_2,start_out2:end_out2) = frame_in(idx_2,start_in2:end_in2);
    
    Y_out(:,:,n_fr) = frame_out;
end

end