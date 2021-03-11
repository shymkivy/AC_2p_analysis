function data = f_s1_load_ca_data(data, ops)
% load average ca traces for alignment
% if multiple voltage files are provided, the ca is assumed to be
% concatenated.
%num_frames = data.num_frames;
%frame_data = data.frame_data;

% ca videos are assumed to be all concatinated into one
if strcmp(ops.ca_processing ,'onacid')
    data.file_cuts_params = cell(1,ops.num_planes);
    data.ave_trace = cell(1,ops.num_planes);
    for n_pl = 1:ops.num_planes
        temp_params = load([ops.file_path_full_ca{n_pl} '_h5cutsinfo']);
        data.file_cuts_params{n_pl} = temp_params.params;
        data.ave_trace{n_pl} = if_normalize(data.file_cuts_params{n_pl}.ave_trace);
    end
elseif strcmp(ops.ca_processing ,'raw_movie')
    Y = double(bigread3([ops.file_dir '\' ops.file_core '.tif']));
    data.ave_trace = cell(1,ops.num_planes);
    trace = squeeze(mean(mean(Y,1),2));
    data.ave_trace{1} = if_normalize(trace);
else
    ops.halo_sub = 1;
    
    clicktrate_data = cell(1,ops.num_planes);
    data.ave_trace = cell(1,ops.num_planes);
    for n_pl = 1:ops.num_planes
        temp_trace = dlmread([ops.file_dir '\' 'TracesRAW_',ops.file_core,'_REDCHAN.csv'],',',0,0)';
        data.ave_trace{n_pl}=if_normalize(temp_trace);
        
        clicktrate_data{n_pl}.cell_trace = dlmread([ops.file_dir '\' 'TracesRAW_',ops.file_core,'.csv'],',',0,0);
        clicktrate_data{n_pl}.cell_halos = dlmread([ops.file_dir '\' 'halosRAW_',ops.file_core,'.csv'],',',0,0);
        
        if ops.halo_sub
            cell_trace = cell_trace - cell_halos;
        end
        
        clicktrate_data{n_pl}.loc = load(['CONTOURS_' ops.file_name, '.mat']);
        
        % fix if length don't match
        if length(data.ave_trace{n_pl}) ~= sum(data.frame_data.num_frames_mpl(:,n_pl))
            ops.errors = [ops.errors; {'Frame numbers from XML and Clicktrace are unequal, Line 164ish'}];
            % fix here
            clicktrate_data{n_pl}.cell_trace = clicktrate_data{n_pl}.cell_trace(:,1:sum(data.frame_data.num_frames_mpl(:,n_pl)));
            data.ave_trace{n_pl} = data.ave_trace{n_pl}(1:sum(data.frame_data.num_frames_mpl(:,n_pl)));
        end
    end
    data.clicktrace_data = clicktrace_data;
    
end

% check if file sizes are weird
for n_pl = 1:ops.num_planes
    num_frames_xml = sum(data.frame_data.num_frames_mpl(:,n_pl));
    num_frames_ca_trace = numel(data.ave_trace{n_pl});
    if num_frames_ca_trace > num_frames_xml
        ops.errors = [ops.errors; {['Fewer XML frames than in movies num_plane ' num2str(n_pl)]}];
        temp_trace = data.ave_trace{n_pl};
        data.ave_trace{n_pl} = temp_trace(1:num_frames_xml);
    elseif num_frames_ca_trace < num_frames_xml
        ops.errors = [ops.errors; {['More XML frames than in movies num_plane ' num2str(n_pl)]}];
%         crop_files = ops.num_volt_in_files;
%         % fix the by cropping at the end of each file
%         extra_XML_frames = num_frames_xml - num_frames_ca_trace;
%         for n_fr = 1:extra_XML_frames
%             temp_inx = rem(n_fr-1,crop_files)+1;
%             data.frame_data.num_frames_mpl(temp_inx,n_pl) = data.frame_data.num_frames_mpl(temp_inx,n_pl) - 1;
%             data.frame_data.num_frames(temp_inx) = data.frame_data.num_frames(temp_inx) - 1;
%             data.frame_data.frame_times_mpl{temp_inx,n_pl}(end) = [];
%             data.frame_data.frame_times{temp_inx}(end) = [];
%         end
    end
end

% combine over multiplane data
data.ave_trace_superpos = f_s1_multiplane_combine(data.ave_trace);

% divide into submovies if was preconcatenated
data.ave_trace_superpos_parts = cell(ops.num_volt_in_files,1);
current_start_frame = 1;
for n_file = 1:ops.num_volt_in_files
    num_frames = data.frame_data_XML.num_frames_raw(n_file);
    data.ave_trace_superpos_parts{n_file} = data.ave_trace_superpos(current_start_frame:(current_start_frame+num_frames-1));
    current_start_frame = current_start_frame + data.frame_data_XML.num_frames_raw(n_file);
end

end


function norm_trace = if_normalize(trace)

base = min(trace);
base_sub = trace - base;
peak = max(base_sub);
norm_trace = base_sub/peak;

end
