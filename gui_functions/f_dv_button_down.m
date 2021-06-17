function f_dv_button_down(app, src)

% get coordinates of mouse click and type of click
if ~isempty(app.DatasetDropDown.Value)
    n_pl = app.mplSpinner.Value;
    est1 = app.ddata.OA_data{n_pl}.est;
    proc1 = app.ddata.OA_data{n_pl}.proc;
    
    
    
    info = get(src);
    coord = round(info.Parent.CurrentPoint(1,1:2));
    indx_current =  sub2ind(proc1.dims', coord(2), coord(1));
    %selection_type = app.UIFigure.SelectionType;
    %app.last_cell_num = app.current_cell_num;
    
    pix_vals = est1.A(indx_current,app.cdata.accepted_cells);
    [temp_val, n_cell] = max(pix_vals);
    if full(temp_val) > 0
  
        contours_accepted = est1.contours(app.cdata.accepted_cells);
        temp_contours = contours_accepted{n_cell};

        if isgraphics(app.gui_plots.plot_current_contour)
            delete(app.gui_plots.plot_current_contour);
        end

        hold(app.UIAxes, 'on');
        app.gui_plots.plot_current_contour = plot(app.UIAxes, temp_contours(:,1), temp_contours(:,2), 'color', [0.75, 0, 0.75], 'LineWidth', 2);
        hold(app.UIAxes, 'off');
        
        app.CellSpinner.Value = n_cell;
        f_dv_update_cell(app);
    end
end

end
