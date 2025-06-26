function f_dv_some_cell_plots(app)

ddata = app.ddata;
cdata = ddata.cdata{1};
snr = ddata.OA_data{1}.proc.SNR2_vals;

min_snr_thresh = prctile(snr, 80);
idx1 = snr > min_snr_thresh;

raw2 = cdata.raw(idx1,:);
S = cdata.S_sm(idx1,:);

[num_cells, num_t] = size(raw2);

start_t = 1000;
len = 10000;
num_to_plot = 10;
plot_t = (1:len)/cdata.volume_period;

cells_idx = randsample(num_cells, num_to_plot);

figure(render='painters'); 
subplot(2,1,1); hold on; axis tight
for n_cell = 1:num_to_plot
    trace1 = raw2(cells_idx(n_cell),start_t:(start_t+len-1));
    trace1 = trace1/max(trace1);
    plot(plot_t, trace1 - 1*n_cell + num_to_plot)
end
xlabel('Time (s)');
ylabel('cells')
title('Raw traces');

subplot(2,1,2); hold on; axis tight
for n_cell = 1:num_to_plot
    trace1 = S(cells_idx(n_cell),start_t:(start_t+len-1));
    trace1 = trace1/max(trace1);
    plot(plot_t, trace1 - 1*n_cell + num_to_plot)
end
xlabel('Time (s)');
ylabel('cells')
title('Deconvolved traces');

accc = ddata.OA_data{1}.proc.comp_accepted;
A  = reshape(sum(ddata.OA_data{1}.est.A(:,accc),2), 256, 256);
figure; imagesc(A'); axis equal tight off;

end