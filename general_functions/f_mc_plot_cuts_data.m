function f_mc_plot_cuts_data(cuts_data, title_tag)

num_planes = numel(cuts_data);
num_it = numel(cuts_data{1}.dsall);
colors1 = parula(num_planes);

[dsall1, dsall1_all, dsall1_all_mf] = f_mc_dsall_proc(cuts_data);

it_disp = zeros(num_it,num_planes);
figure; hold on;
for n_pl = 1:num_planes
    for n_it = 1:num_it
        it_disp(n_it, n_pl) = mean(sqrt(sum(cuts_data{n_pl}.dsall{n_it}.^2,2)));
    end
end
plot(it_disp, '-o'); xlabel('iteration'); ylabel('mean disp');
title([title_tag, ' mean fix'], 'interpreter', 'none');

for n_pl = 1:num_planes
    figure;
    sp_all = cell(num_it+1,1);
    for n_it = 1:num_it
        sp_all{n_it} = subplot(num_it+1, 1, n_it);
        plot(cuts_data{n_pl}.dsall{n_it});
        axis tight;
    end
    sp_all{num_it + 1} = subplot(num_it+1, 1, num_it+1);
    plot(sum(cat(3,cuts_data{n_pl}.dsall{:}),3));
    axis tight;
    ylabel('total fix');
    linkaxes([sp_all{:}], 'x');  axis tight;
    sgtitle(sprintf('Fix per iteration; pl%d; %s', n_pl, title_tag), 'interpreter', 'none')
end

figure;
sp1 = subplot(2,1,1); hold on;
for n_pl = 1:num_planes
    plot(dsall1{n_pl}(:,1), 'color', colors1(n_pl,:));
end
plot(dsall1_all(:,1), 'k');
plot(dsall1_all_mf(:,1), 'g');
ylabel('y motion');
title(sprintf('%s moco',title_tag), 'interpreter', 'none');
sp2 = subplot(2,1,2); hold on;
for n_pl = 1:num_planes
    plot(dsall1{n_pl}(:,2), 'color', colors1(n_pl,:));
end
plot(dsall1_all(:,2), 'k');
plot(dsall1_all_mf(:,2), 'g');
ylabel('x motion');
linkaxes([sp1 sp2], 'x');  axis tight;




end