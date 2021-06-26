function firing_rate = f_dv_get_firing_rate(cdata)

smooth_SD = 0; % 110 is better?

%%
volume_period = cdata.volume_period;

firing_rate = cdata.S;

%%
active_cells = sum(firing_rate,2) ~= 0;
firing_rate(~active_cells,:) = [];

num_cells = size(firing_rate,1);

firing_rate = firing_rate(randperm(num_cells),:);

firing_rate = f_smooth_gauss(firing_rate, smooth_SD/volume_period);

end