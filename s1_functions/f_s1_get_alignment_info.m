function data = f_s1_get_alignment_info(data, ops)

% ask if voltage traces need to be aligned
if ~isfield(data, 'alignment')
    figure; hold on;
    plot(1:data.num_samp_volt,data.volt_data_all(:,ops.align_to_channel));
    plot((1:data.frame_data.num_frames_all)*data.frame_data.frame_period_ave,data.ave_trace_superpos);
    title(sprintf('Do these need to be aligned? [Y/N];\n volt trig mode: %s', data.frame_data_XML.trigger_mode{1}))
    alignment.need_alignment = ask_yes_no_fig();
    close;
    
    if alignment.need_alignment
        if sum(strcmp(ops.alignment_method, {'xcorr','peak_onsets','manual', 'peak_onsets_scale_only', 'peak_onsets_shift_only'}))
            alignment_method = ops.alignment_method;
        else
            alignment_method = 'xcorr'; % default
        end

        alignment.align_method = cell(ops.num_volt_in_files,1);
        alignment.scaling_factor = zeros(ops.num_volt_in_files,1);
        alignment.shift = zeros(ops.num_volt_in_files,1);
        alignment.alignment_good = zeros(ops.num_volt_in_files,1);

        for n_file = 1:ops.num_volt_in_files
            volt_data = data.volt_data_all(:,ops.align_to_channel);
            ca_traces = data.ave_trace_superpos_parts{n_file};
            frame_times = data.frame_data_XML.frame_times_raw{n_file};
            temp_method = alignment_method;

            temp_finish_alignment = 0;
            align_now = 1;
            figure;
            while ~temp_finish_alignment

                if align_now
                    if strcmp(temp_method, 'xcorr')
                        shift = f_s1_align_traces_xcor(ca_traces, frame_times, volt_data);
                        scaling_factor = 1;
                    elseif strcmp(temp_method, 'peak_onsets')
                        [shift, scaling_factor] = f_s1_align_traces_peak_onesets(ca_traces, frame_times, volt_data);
                    elseif strcmp(temp_method, 'peak_onsets_shift_only')
                        [shift, ~] = f_s1_align_traces_peak_onesets(ca_traces, frame_times, volt_data);
                        scaling_factor = 1;
                    elseif strcmp(temp_method, 'peak_onsets_scale_only')
                        [~, scaling_factor] = f_s1_align_traces_peak_onesets(ca_traces, frame_times, volt_data);
                        shift = 0;
                    elseif strcmp(temp_method, 'manual')
                        [shift, scaling_factor] = f_s1_align_traces_manual(ca_traces, frame_times, volt_data);
                    else
                        error('Select proper alignment method')
                    end
                end

                % try aligning
                shifted_scaled_volt_data = align_volt_by_scale_shift2(volt_data, scaling_factor, shift);

                %verify the alignment

                clf;
                plot(1:(size(shifted_scaled_volt_data,1)),shifted_scaled_volt_data);
                hold on;
                plot(frame_times, ca_traces);
                title(sprintf('Aligned traces with %s; Is it good? [Y/N]\n shift = %.2fms; scale = %.7f', temp_method,shift,scaling_factor), 'Interpreter', 'none');
                legend('DAQ voltage trace', 'Alignment Channel');
                fprintf('Are the alignments good? [Y/N](click on fig to answer)\n');
                % if alignment is ok, keep manual alignment as 1
                alignment_good = ask_yes_no_fig();

                if alignment_good
                    temp_finish_alignment = 1;
                else
                    ans_temp = if_get_input('Enter 1 for xcorr; 2 for peak_onsets; 3 for manual alignment; 4 for manual shift; 0 to not align: ', [0,1,2,3,4]);
                    if ans_temp == 1
                        temp_method = 'xcorr';
                        align_now = 1;
                    elseif ans_temp == 2
                        temp_method = 'peak_onsets';
                        align_now = 1;
                    elseif ans_temp == 3
                        temp_method = 'manual';
                        align_now = 1;
                    elseif ans_temp == 4
                        fine_shift_adj = if_get_input('Enter number of frames to shift the voltage to the right: ', 'numeric');
                        shift = shift + round(fine_shift_adj);
                        align_now = 0;
                    elseif ~ans_temp
                        temp_finish_alignment = 1;
                        shift = 0;
                        scaling_factor = 1;
                        temp_method = 0;
                    end
                end
            end 
            close;
            alignment.scaling_factor(n_file) = scaling_factor;
            alignment.shift(n_file) = shift;
            alignment.align_method{n_file} = temp_method;
            alignment.alignment_good(n_file) = alignment_good;
        end

        % check if some of the traces were not manually aligned, and
        % compute thaverage alignment parameters 
        if sum(alignment.alignment_good)
            for n_file = 1:ops.num_volt_in_files
                if ~alignment.alignment_good(n_file)
                    figure;
                    plot(1:(size(shifted_scaled_volt_data,1)),shifted_scaled_volt_data);
                    hold on;
                    plot(frame_times, ca_traces);
                    title(['Alignment ' num2str(n_file) ' was skipped, do you want to apply average alignment from other traces? [Y/N]']);
                    legend('DAQ voltage trace', 'Alignment Channel');

                    % if alignment is ok, keep manual alignment as 1
                    apply_mean = ask_yes_no_fig();
                    close;
                    if apply_mean
                        alignment.scaling_factor(n_file) = mean(alignment.scaling_factor(alignment.alignment_good));
                        alignment.shift(n_file) = mean(alignment.shift(alignment.shift));
                        alignment.align_method{n_file} = 'average_of_others';
                    end
                end
            end
        end
        
    else
        alignment.shift = 0;
        alignment.scaling_factor = 1;
    end
    data.alignment = alignment;
    if exist(ops.file_save_path_full_processing_params, 'file')
        save(ops.file_save_path_full_processing_params, 'alignment', '-append')
    else
        save(ops.file_save_path_full_processing_params, 'alignment')
    end
    fprintf('Saved alignment data\n');
end

end

function inp_answer = if_get_input(question, outputs)

got_answer = 0;
while ~got_answer
    inp_answer = input(question);
    if strcmp(outputs,'numeric')
        if isnumeric(inp_answer)
            got_answer = 1;
        else
            fprintf(['Choose a valid option [' outputs ']:\n']);
        end
    else
        if sum(inp_answer == outputs)
            got_answer = 1;
        else
            fprintf(['Choose a valid option [' num2str(outputs) ']:\n']);
        end
    end
end

end
