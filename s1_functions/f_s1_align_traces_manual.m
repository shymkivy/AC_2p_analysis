function [shift, scaling_factor] = f_s1_align_traces_manual(ca_traces, frame_times, volt_data)
v_times = 1:size(volt_data,1);

temp_finish_alignment = 0;
while temp_finish_alignment == 0

    % plot data to find alignment parameters
    figure;
    dcm_obj = datacursormode;
    set(dcm_obj,'UpdateFcn', @NewCallback_YS)
    plot(v_times,volt_data);
    hold on;
    plot(frame_times, ca_traces);
    title('Traces before correction, enter coordinates in prompt');
    legend('DAQ voltage trace', 'Alignment channel');

    % ask how many points to use for alignment
    num_synch_points = input('How many alignment points to use? (0 to cancel)');


    if num_synch_points > 0

        x = zeros(2,num_synch_points);
        for ii = 1:num_synch_points
            x(1,ii) = input(sprintf('Input synch point %d of Ca trace (top):\n', ii));
            x(2,ii) = input(sprintf('Input synch point %d of Voltage trace (top):\n', ii));
        end

        % extract the alignment information
        if size(x,2) > 1
            scaling_factor = (x(1,2)-x(1,1))/(x(2,2)-x(2,1));
        else
            scaling_factor = 1;
            warning('Cannot get scale factor with just one point');
        end
        shift = round(x(1,1) - (x(2,1)*scaling_factor));

        close;

        % try aligning
        shifted_scaled_dat_proc = align_volt_by_scale_shift2(volt_data, scaling_factor, shift);

        %verify the alignment
        figure;
        plot(1:numel(shifted_scaled_dat_proc),shifted_scaled_dat_proc);
        hold on;
        plot(frame_times, ca_traces);
        title('Aligned traces; Is it good? [Y/N]');
        legend('DAQ voltage trace', 'Alignment Channel');

        % if alignment is ok, keep manual alignment as 1
        temp_finish_alignment = ask_yes_no_fig();

    else
        temp_finish_alignment = 1;
        shift = 0;
        scaling_factor = 1;
    end
    close;
end

end