function FactorAxes = f_subplot_tight(m,n,idx,data_t,data_x, varargin)
%VIZ_KTENSOR Visualize a ktensor
%
% Xticks = cell array (one per mode)
% Xticklabels = cell array (one per mode)
% Xlims = cell array (one per mode)
% Ylims = cell array (one per mode)
% Hrelspace = Relative Horizontal Space Between Plots (.1)
% HspaceLeft = Absolute Space on Left of Each Plot (0)
% HspaceRight = Absolute Space on right of Each Plot (0)
% Vrelspace = Relative Hoizontal Space Between Plots (.1)
% VspaceTop = Absolute Space for Title (.1)
% VspaceBottom = Absolute Space for XTick Marks (.1)

%%
% horizontal
[n_idx, m_idx] = ind2sub([n, m],idx);
% vertical
%[m_idx, n_idx] = ind2sub([m, n],idx);


% parse optional inputs

params = inputParser;
% Figure 
params.addParameter('Figure', []);
% Spacing
params.addParameter('Relmodespace', ones(n,1)); % Horizontal space for each mode
params.addParameter('Hspace',0.01); % Horizontal space between axes
params.addParameter('Hspaceright',0.025); % Horizontal space on left
params.addParameter('Hspaceleft',0.05); % Horizontal space on right
params.addParameter('Vspace',0.01); % Vertical space between axes
params.addParameter('Vspacetop',0.05); % Vertical space at top
params.addParameter('Vspacebottom',0.05); % Vertical space at bottom
% Titles
params.addParameter('Modetitles', []);
params.addParameter('Factortitle', 'none'); % Default is 'none'. Options are 'weight' or 'number'
% Plots
params.addParameter('Plottype', repmat({'line'}, [n 1]));
params.addParameter('Plotsize', -1 * ones(n,1)); % Used for scatter dot size or plot linewidth
params.addParameter('Plotcolors', cell(n,1));
params.addParameter('Sameylims', true(m,1));


params.parse(varargin{:});
res = params.Results;

%% Create figure
% if isempty(res.Figure)
%     figure;
% else
%     figure(res.Figure);
%     clf;
% end

%% Create axes
Vplotspace = 1 - res.Vspacetop - res.Vspacebottom - (m - 1) * res.Vspace;
height = Vplotspace / m;

Hplotspace = 1 - res.Hspaceleft - res.Hspaceright - (n - 1) * res.Hspace;
width = (res.Relmodespace ./ sum(res.Relmodespace)) .* Hplotspace;

% Global axis
GlobalAxis = axes('Position',[0 0 1 1]); % Global Axes
axis off;



% Factor axes
FactorAxes = gobjects(1,1); % Factor Axes

xpos = res.Hspaceleft + (n_idx-1) * res.Hspace + sum(width(1:n_idx-1));
ypos = 1 - res.Vspacetop - height - (m_idx-1) * (height + res.Vspace);
FactorAxes = axes('Position',[xpos ypos width(n_idx) height]);
%set(FactorAxes,'FontSize',14);



% %% Plot each factor
% h = gobjects(1,1);
% if res.Plotsize(n_idx) == -1
%     lw = 1;
%     ss = 10;
% else
%     lw = res.Plotsize(n_idx);
%     ss = res.Plotsize(n_idx);
% end
% 
% if isempty(res.Plotcolors{n_idx})
%     cc = [0 0 1];
% else
%     cc = res.Plotcolors{n_idx};
% end
% 
% 
% xl = [min(data_t) max(data_t)]; %data_t
% yl = [min(data_x) max(data_x)];
%        
% xx = data_t;
% yy = data_x;
% hold(FactorAxes, 'off');
% 
% switch res.Plottype{n_idx}
%     case 'line'
%         hh = plot(FactorAxes, xx, yy, 'Linewidth', lw, 'Color', cc);
%     case 'scatter'
%         hh = scatter(FactorAxes, xx, yy, ss, cc, 'filled');
%     case 'bar'
%         hh = bar(FactorAxes, xx, yy, 'EdgeColor', cc, 'FaceColor', cc);
% end
% 
% xlim(FactorAxes,xl);
% if res.Sameylims(m_idx)
%     ylim(FactorAxes,yl);
% else
%     tmpyl = yl;
%     ylim(FactorAxes,tmpyl);
% end
% set(FactorAxes,'Ytick',[]);
% if m_idx < m
%     set(FactorAxes,'XtickLabel',{});
% end            
% 
% hold(FactorAxes, 'on');
% plot(FactorAxes, xl, [0 0], 'k:', 'Linewidth', 1.5);
% 
% h = hh;
% set(FactorAxes,'FontSize',14)


% %% Title for each mode
% htitle = gobjects(n,1);
% if ( isscalar(res.Modetitles) && islogical(res.Modetitles) && (res.Modetitles == false) )
%     ModeTitles = 'none';
% else
%     if isempty(res.Modetitles)
%         ModeTitles = cell(n,1);
%         for i = 1:n
%             ModeTitles{i} = sprintf('Mode %d',i);
%         end
%     else
%         ModeTitles = res.Modetitles;
%     end
%     
%     axes(GlobalAxis);
%     for n_idx = 1:n
%         xpos = res.Hspaceleft + (n_idx-1) * res.Hspace + sum(width(1:n_idx-1)) + 0.5 * width(n_idx);
%         %xpos = res.Hspaceleft + (k-1) * (width + res.Hspace) + 0.5 * width;
%         ypos = 1 - res.Vspacetop;
%         htitle(n_idx) = text(xpos,ypos,ModeTitles{n_idx},'VerticalAlignment','Bottom','HorizontalAlignment','Center');
%         set(htitle(n_idx),'FontSize',16)
%         set(htitle(n_idx),'FontWeight','bold')
%     end
% end
% 
% %% Print factor titles
% hftitle = gobjects(m,1);
% if ~strcmpi(res.Factortitle,'none')
%     axes(GlobalAxis);
%     rellambda = abs (K.lambda / K.lambda(1));
%     for m_idx = 1:m
%         xpos = 0.9 * res.Hspaceleft;
%         ypos = 1 - res.Vspacetop - 0.5 * height - (m_idx-1) * (height + res.Vspace);
%         %ypos = 1 - res.Vspacetop - 0.5 * height - (j-1) * (1 + res.Vrelspace) * height;
%         if strcmpi(res.Factortitle,'weight')          
%             txt = sprintf('%3.2f', rellambda(m_idx));
%         else
%             txt = sprintf('%d', m_idx);
%         end
%         hftitle(m_idx) = text(xpos,ypos,txt,'VerticalAlignment','Middle','HorizontalAlignment','Right');
%         set(hftitle(m_idx),'FontSize',14)
%     end
% end
% %% Save stuff to return
% info.height = height;
% info.width = width;
% info.ModeTitles = ModeTitles;
% info.GlobalAxis = GlobalAxis;
% info.FactorAxes = FactorAxes;
% info.htitle = htitle;
% info.hftitle = hftitle;
% info.h = h;
% 
% axes(GlobalAxis)
