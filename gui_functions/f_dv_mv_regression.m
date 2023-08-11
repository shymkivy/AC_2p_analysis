function f_dv_mv_regression(app)

disp("don't work")

n_pl = 1;

ddata = app.ddata;

stim_regress_duration = 1000;

vol_period = ddata.proc_data{1}.frame_data.volume_period;

regress_frames = round(stim_regress_duration/vol_period);


tn_all = f_dv_get_trial_number(app);
tt = app.ops.context_types_all(tn_all)';

tt = 1:10;
num_tt = numel(tt);

cdata = f_dv_get_cdata(app);
firing_rate = cat(1,cdata.S_sm);

[num_cells, num_t] = size(firing_rate);




stim_trace = zeros(num_t,num_tt, regress_frames);
for n_tt = 1:num_tt
    idx_stim = app.ddata.stim_frame_index{n_pl}(logical(sum(ddata.trial_types{1} == tt(n_tt),2)));
    for n_stim = 1:numel(idx_stim)
        for n_reg = 1:regress_frames
            stim_trace(idx_stim(n_stim)+regress_frames-1, n_tt, n_reg) = 1;
        end
    end
end

stim_trace2 = reshape(stim_trace, num_t, []);


figure;
plot(stim_trace)


[beta,Sigma,E,CovB,logL] = mvregress(stim_trace2, firing_rate(1, :)');


figure; plot(beta)


x = stim_trace2*beta;

figure; hold on;
plot(x)
plot(firing_rate(1, :))



end