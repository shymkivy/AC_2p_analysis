function norm_trace = f_normalize_YS(trace)

base = min(trace);
base_sub = trace - base;
peak = max(base_sub);
norm_trace = base_sub/peak;

end