function f_dv_plot_cont_dev_red(app)

bar_sm_ylim = [0 .85];
bar_dist_ylim = [0 7];


[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

resp_cell_type = app.ResposivecellstypeDropDown.Value;

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

num_dsets = size(data,1);

tn_all2 = [18 19 20; 28 29 30];

data_all = cell(num_dsets,1);
for n_dset = 1:num_dsets
    data1 = cell(2,1);
    for n_fl = 1:2
        stats1 = cat(1,data(n_dset,:).stats{n_pl});
        
        if strcmpi(app.ResponsivecellsselectDropDown.Value, 'All')
            resp_cell_sel = 'All';
        else
            resp_cell_sel = 'Resp marg';
        end

        [~, resp_vals] = f_dv_get_resp_vals_cells(app, stats1, tn_all2(n_fl,:), [], resp_cell_sel);
        
        data1{n_fl} = cat(2,resp_vals{:});
    end
    data_all{n_dset} = cat(1,data1{:});
end

data_all2 = cat(1,data_all{:});

data_all3 = data_all2;

max_pt = max(data_all2(:));
min_pt = round(max_pt*.1,2);

title_tag2 = sprintf('%s; resp %s; %s', resp_cell_type, resp_cell_sel, title_tag);

f1 = figure;
plot3(data_all2(:,1), data_all2(:,2), data_all2(:,3), '.k', 'MarkerSize', 1)
xlabel('cont'); ylabel('red'); zlabel('dev');
f1.Children.XDir = 'reverse';
f1.Children.YDir = 'reverse';
f1.Children.XLim = [-min_pt max_pt+min_pt];
f1.Children.YLim = [-min_pt max_pt+min_pt];
f1.Children.ZLim = [-min_pt max_pt+min_pt];
grid on
title(sprintf('%s', title_tag2), 'interpreter', 'none')


plor_pairs = [1, 3; 2, 3; 2, 1];
title_ctx = {'cont', 'red', 'dev'};

num_pl = size(plor_pairs,1);

cos_si_all = zeros(num_pl,1);
dist_all = zeros(num_pl,1);
bar_leg = cell(3,1);
for n_pl = 1:num_pl
    pp2 = plor_pairs(n_pl,:);
    
    data_all3(or(isnan(data_all2(:,pp2(1))), isnan(data_all2(:,pp2(2)))),:) = [];
    cos_sim1 = 1 - pdist2(data_all3(:,pp2(1))', data_all3(:,pp2(2))', 'cosine');
    dist1 = pdist2(data_all3(:,pp2(1))', data_all3(:,pp2(2))', 'euclidean');
    
    cos_si_all(n_pl) = cos_sim1;
    dist_all(n_pl) = dist1;
    bar_leg{n_pl} = sprintf('%s vs %s', title_ctx{pp2(1)}, title_ctx{pp2(2)});
    
    
    f2 = figure;
    plot(data_all2(:,pp2(1)), data_all2(:,pp2(2)), '.k', 'MarkerSize', 1);
    xlabel(title_ctx{pp2(1)}); ylabel(title_ctx{pp2(2)});
    f2.Children.XLim = [-min_pt max_pt+min_pt];
    f2.Children.YLim = [-min_pt max_pt+min_pt];
    title(sprintf('%s; %s, cos SI = %.2f', title_tag2, bar_leg{n_pl}, cos_sim1), 'interpreter', 'none')
end


figure;
subplot(1,2,1);
bar(categorical(categorical(bar_leg)), cos_si_all)
ylim(bar_sm_ylim);
ylabel('Cosine similarity');
title('Cosine similarity');
subplot(1,2,2);
bar(categorical(categorical(bar_leg)), dist_all);
ylim(bar_dist_ylim);
title('Euclidean distance');
ylabel('Euclidean distance');
sgtitle(sprintf('%s', title_tag2));
% 
% f3 = figure;
% plot(data_all2(:,2), data_all2(:,3), '.k', 'MarkerSize', 1);
% xlabel('red'); ylabel('dev');
% f3.Children.XLim = [-min_pt max_pt+min_pt];
% f3.Children.YLim = [-min_pt max_pt+min_pt];
% title(sprintf('%s; red vs dev', title_tag2), 'interpreter', 'none')
% 
% f4 = figure;
% plot(data_all2(:,2), data_all2(:,1), '.k', 'MarkerSize', 1);
% xlabel('red'); ylabel('cont');
% f4.Children.XLim = [-min_pt max_pt+min_pt];
% f4.Children.YLim = [-min_pt max_pt+min_pt];
% title(sprintf('%s; red vs cont', title_tag2), 'interpreter', 'none')

end