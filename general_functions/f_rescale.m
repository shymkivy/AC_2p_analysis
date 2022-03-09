function trace_out = f_rescale(trace_in)

base = min(trace_in);
base_sub = trace - base;
peak = max(base_sub);
trace_out = base_sub/peak;

end