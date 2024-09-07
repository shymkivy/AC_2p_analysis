function f_plot_bar_points(data, labels, colors_big, colors_small, error_type, do_violin)

if ~exist('error_type', 'var') || isempty(error_type)
    error_type = 'sem';
end

if ~exist('do_violin', 'var') || isempty(do_violin)
    do_violin = 1;
end

[num_sbar, num_bbar] = size(data);

means1 = zeros(num_sbar, num_bbar);
sems1 = zeros(num_sbar, num_bbar);
for n_bb = 1:num_bbar
    num_pts = size(data{1,n_bb},1);
    for n_sb = 1:num_sbar
        temp_data = data{n_sb,n_bb};
        means1(n_sb,n_bb) = mean(temp_data,1);
        if strcmpi(error_type, 'sem')
            sems1(n_sb,n_bb) = std(temp_data,[],1)/sqrt(num_pts-1);
        elseif strcmpi(error_type, 'std')
            sems1(n_sb,n_bb) = std(temp_data,[],1);
        end
    end
end

figure; hold on;
bar(categorical(labels(:,1), labels(:,1)), zeros(num_bbar,1));

for n_bb = 1:num_bbar
    if num_sbar > 1
        locs1 = n_bb + ((1:4)-2.5)/5.5;
        siz_spread = 1/15;
        for n_sb = 1:num_sbar
            temp_data = data{n_sb,n_bb};
            if do_violin
                [x1, y1] = if_violin_points(temp_data);
            else
                x1 = (rand(numel(temp_data),1)-0.5);
                y1 = temp_data;
            end
            plot(x1*siz_spread+locs1(n_sb), y1, '.', color=colors_small{n_sb})
        end
    else
        locs1 = n_bb;
        siz_spread = 1/2;
        temp_data = data{1,n_bb};
        if do_violin
            [x1, y1] = if_violin_points(temp_data);
        else
            x1 = (rand(numel(temp_data),1)-0.5);
            y1 = temp_data;
        end
        plot(x1*siz_spread+locs1(n_sb), y1, '.', color=colors_big{n_bb}, Markersize=1);
    end
    b1 = bar(n_bb, means1(:,n_bb));
    if num_sbar > 1
        errorbar(locs1, means1(:,n_bb), sems1(:,n_bb), '.k');
        for n_sb = 1:num_sbar
            b1(n_sb).FaceColor = colors_small{n_sb};
            b1(n_sb).FaceAlpha = 0.4;
        end
    else
        errorbar(n_bb, means1(:,n_bb), sems1(:,n_bb), '.k');
        b1.FaceColor = colors_big{n_bb};
        b1.FaceAlpha = 0.4;
    end
end


end

function [x1, y1] = if_violin_points(y0)

y1 = sort(y0);
x0 = (rand(numel(y0),1)-0.5);
x1 = x0;
num_bins = 100;
counts = histcounts(y0,num_bins);
countsn = counts/max(counts);
idx_start = 1;
for n_bin = 1:num_bins
    idx_end = idx_start + counts(n_bin) - 1;
    x1(idx_start:idx_end) = x0(idx_start:idx_end) .* countsn(n_bin);
    idx_start = idx_end+1;
end
end