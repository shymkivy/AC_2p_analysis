function f_dv_button_down(app, src)

% get coordinates of mouse click and type of click
if ~isempty(app.DatasetDropDown.Value)
    n_pl = app.mplSpinner.Value;
    est1 = app.ddata.OA_data{n_pl}.est;
    proc1 = app.ddata.OA_data{n_pl}.proc;
    accepted_cells = app.cdata{n_pl}.accepted_cells;
    
    info = get(src);
    coord = round(info.Parent.CurrentPoint(1,1:2));
    indx_current =  sub2ind(proc1.dims', coord(2), coord(1));
    %selection_type = app.UIFigure.SelectionType;
    %app.last_cell_num = app.current_cell_num;
    
    A_acc = est1.A(:,accepted_cells);
    acc_idx = find(accepted_cells);
    
    pix_vals = A_acc(indx_current,:);
    [temp_val, n_cell] = max(pix_vals);
    if full(temp_val) > 0
  
        contours_accepted = est1.contours(accepted_cells);
        temp_contours = contours_accepted{n_cell};

        if isgraphics(app.gui_plots.plot_current_contour)
            delete(app.gui_plots.plot_current_contour);
        end

        hold(app.UIAxes, 'on');
        app.gui_plots.plot_current_contour = plot(app.UIAxes, temp_contours(:,1), temp_contours(:,2), 'color', [0.75, 0, 0.75], 'LineWidth', 2);
        hold(app.UIAxes, 'off');
        
        app.CellSpinner.Value = n_cell;
        f_dv_update_cell(app);
        
        A_im = reshape(A_acc(:,n_cell), proc1.dims');
        app.gui_plots.image_roi.CData = A_im;
        title(app.UIAxes_roi, sprintf('Cell %d', n_cell));
        
        y_lim = find(sum(A_im,2)>0);
        y_coord = [y_lim(1)-1, y_lim(end)+1];
        y_size = y_coord(2) - y_coord(1);
        x_lim = find(sum(A_im,1)>0);
        x_coord = [x_lim(1)-1, x_lim(end)+1];
        x_size = x_coord(2) - x_coord(1);
        xy_diff = floor(abs(y_size - x_size)/2);
        if y_size > x_size
            xlim(app.UIAxes_roi, x_coord + [-xy_diff, xy_diff]);
            ylim(app.UIAxes_roi, y_coord);
        elseif y_size < x_size
            xlim(app.UIAxes_roi, x_coord);
            ylim(app.UIAxes_roi, y_coord + [-xy_diff, xy_diff]);
        else
            xlim(app.UIAxes_roi, x_coord);
            ylim(app.UIAxes_roi, y_coord);
        end
        
        app.SNRcaimanEditField.Value = est1.SNR_comp(acc_idx(n_cell));
        app.SNR2EditField.Value = proc1.SNR2_vals(acc_idx(n_cell));
        app.CNNprobEditField.Value = est1.cnn_preds(acc_idx(n_cell));
        app.RvalueEditField.Value = est1.r_values(acc_idx(n_cell));
    end
end

end
