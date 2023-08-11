function tn_out = f_dv_proc_tn(tn_in, mmn_freq)

tn_out = tn_in;

idx1 = and(tn_in > 30, tn_in < 50);
tn_out(idx1) = tn_out(idx1) - 40 + mmn_freq(2);

idx1 = and(tn_in > 50, tn_in < 70);
tn_out(idx1) = tn_out(idx1) - 60 + mmn_freq(1);

end