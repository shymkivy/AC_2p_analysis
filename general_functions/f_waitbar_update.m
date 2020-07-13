function f_waitbar_update(wb, fraction_comp)

if wb.new_fig
    waitbar(fraction_comp,wb.h);
else
    wb.handlew.Value = fraction_comp;
end

end