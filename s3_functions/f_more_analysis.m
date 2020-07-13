function f_more_analysis(data, ops)
    for n_cond = 1:numel(ops.conditions_to_analyze)
        cond_name = ops.conditions{ops.conditions_to_analyze(n_cond)};

        temp_traces_all = data.(cond_name).vec_data(:,ops.win.no_base_window,:);
        dataset_index = data.(cond_name).num_dataset;
        num_datasets = max(dataset_index);
        
        context_resp_cells = data.(cond_name).context_resp_cells;
        resp_cells = data.(cond_name).resp_cells;
        trial_type = data.(cond_name).trial_type;
        MMN_freq = data.(cond_name).MMN_freq;
        
        context_mmn = data.(cond_name).context_mmn;
        context_mmn_dset = zeros(num_datasets,6);
        trial_type_dset = zeros(num_datasets,size(trial_type,2));
        for n_dset = 1:num_datasets
            context_mmn_dset(n_dset,:) = context_mmn(find(dataset_index==n_dset,1),:);
            trial_type_dset(n_dset,:) = trial_type(find(dataset_index==n_dset,1),:);
        end
        num_cells = data.(cond_name).num_cells;

            

        % extract cells responsive to cont, red, dev, and flip c r d
        context_resp_cells2 = zeros(num_cells,6);
        for n_cell = 1:num_cells
            temp_resp_cells = resp_cells(n_cell,:);
            temp_ctx_index = reshape(context_mmn(n_cell,:,:), 1, []);
            context_resp_cells2(n_cell,:) = temp_resp_cells(temp_ctx_index);
        end
        
        
        % extract percentages of cont/red/dev cells per dataset
        dset_num_cells = zeros(num_datasets,1);
        dset_num_resp_cells = zeros(num_datasets,6);
        dset_num_ctx_intersect = zeros(num_datasets,3);
        dset_num_ctx_trials = zeros(num_datasets,6);
        
        for n_dset = 1:num_datasets
            temp_dset_index = dataset_index == n_dset;
            dset_num_cells(n_dset) = sum(temp_dset_index);
            temp_resp_cells = context_resp_cells2(temp_dset_index,:);
            dset_num_resp_cells(n_dset,:) = sum(temp_resp_cells);
            dset_num_ctx_intersect(n_dset,:) = sum((temp_resp_cells(:,1:3)+temp_resp_cells(:,4:6))==2,1);
            
            for ii = 1:6
                dset_num_ctx_trials(n_dset, ii) = sum(trial_type_dset(n_dset,:)==ops.context_types_all(context_mmn_dset(n_dset,ii)));
            end
        end
        
        dset_means = mean(dset_num_resp_cells./dset_num_cells)*100;
        dset_std = std(dset_num_resp_cells./dset_num_cells)*100;
        dset_intersect_means = mean(dset_num_ctx_intersect./dset_num_cells)*100;
        dset_intersect_std = std(dset_num_ctx_intersect./dset_num_cells)*100;
        
%         x_labels = categorical({'cont1 cont2 intrsc','red1 red2 intrsc','dev1 dev2 intrsc'});
%         data_x = [dset_means(1) dset_means(4) dset_intersect_means(1);...
%                   dset_means(2) dset_means(5) dset_intersect_means(2);...
%                   dset_means(3) dset_means(6) dset_intersect_means(3)];
%         std_x = [dset_std(1) dset_std(4) dset_intersect_std(1);...
%                  dset_std(2) dset_std(5) dset_intersect_std(2);...
%                  dset_std(3) dset_std(6) dset_intersect_std(3)];
%         
%              
%              
%         figure;
%         hBar = bar(x_labels, data_x);
% %         X=cell2mat(get(hBar,'XData'));
% %         hold on;
% %         errorbar(X, data_x,...
% %              data_x - std_x, data_x + std_x);
%         ylabel('Percent of responsive cells');
%         title(sprintf('%s Responsive cells distribution',cond_name));

%         figure;
%         bar(categorical({'cont1 cont2 intrsc','red1 red2 intrsc','dev1 dev2 intrsc'}),...
%             [sum(context_resp_cells2(:,[1 4])) sum(sum(context_resp_cells2(:,[1 4]),2) == 2);...
%             sum(context_resp_cells2(:,[2 5])) sum(sum(context_resp_cells2(:,[2 5]),2) == 2);...
%             sum(context_resp_cells2(:,[3 6])) sum(sum(context_resp_cells2(:,[3 6]),2) == 2);...
%             ]);
%         title([cond_name ' responsive cell counts']);

%         figure;
%         bar(categorical({'cont-dev1','cont-dev2'}),...
%             [sum(sum(context_resp_cells2(:,[1 3]),2) == 2), sum(sum(context_resp_cells2(:,[4 6]),2) == 2)]);
%         title([cond_name ' cells shared between context']);
%         
        
        temp_resp_cells = dataset_index.*context_resp_cells(:,1);
        dataset_resp_cell_counts = zeros(max(dataset_index),1);
        for ii = 1:max(dataset_index)
            dataset_resp_cell_counts(ii) = sum(temp_resp_cells == ii);
        end
        figure;
        bar(dataset_resp_cell_counts);
        title('responding cells per dataset');
        
        % plot for individual datasets the traces
        for n_file = 1:max(dataset_index)
            intr_cells = logical((dataset_index == n_file).*sum(context_resp_cells2(:,1:3),2));
            intr_cells_index = find(intr_cells);
            
            sort_x = unique([find(context_resp_cells2(dataset_index == n_file,3)==1);...
                             find(context_resp_cells2(dataset_index == n_file,1)==1);...
                             find(context_resp_cells2(dataset_index == n_file,2)==1)],'stable');
            
                         
            
            if numel(intr_cells_index)> 0
            
                temp_resp_cells = resp_cells(dataset_index == n_file,:); 

        %         figure;
        %         imagesc(temp_resp_cells)
        %         

                temp_trial_type = trial_type(intr_cells_index(1),:);
                temp_MMN_freq = MMN_freq(intr_cells_index(1),:);

                temp_context_mmn = [temp_MMN_freq(2), 200 + ops.redundent_to_analyze, 170; ...
                                    temp_MMN_freq(1), 100 + ops.redundent_to_analyze, 270]';
                %clear temp_MMN_freq;





                % plot sums of cell acivity
                temp_traces_data = reshape(mean(temp_traces_all(intr_cells,:,:)), 1, []);
                y_max_temp = max(temp_traces_data);
                y_min_temp = min(temp_traces_data);

                figure;
                plot(temp_traces_data);
                hold on;
                for n_cntxt = 1:3
                    temp_ctx_trials = find(logical((temp_trial_type == temp_context_mmn(n_cntxt,1))));
                    temp_ctx_frames = (temp_ctx_trials * ops.win.no_base_window_size) - ops.win.no_base_window_size + 1;
                    for n_frame = 1:numel(temp_ctx_frames)
                        line([temp_ctx_frames(n_frame) temp_ctx_frames(n_frame)], [y_min_temp y_max_temp], 'color', ops.context_colors{n_cntxt});
                    end
                end
                axis tight;
                title([cond_name ' context trials start, dataset ' num2str(n_file)]);
                %clear temp_ctx_trials temp_ctx_frames y_max_temp y_min_temp temp_traces_data;



                % plot each cell dev trial by trial responses
                intr_cells2 = logical((dataset_index == n_file).*context_resp_cells(:,1));
                intr_cells_index2 = find(intr_cells2);



                for n_cntxt = 1:3
                    temp_ctx_trials = find(logical(temp_trial_type == temp_context_mmn(n_cntxt,1)));
                    temp_ctx_traces = temp_traces_all(intr_cells2,:,temp_ctx_trials);

                    temp_ctx_traces_flat_trials = reshape(temp_ctx_traces, numel(intr_cells_index2), []);
                    temp_ctx_traces_flat_cells = reshape(permute(temp_ctx_traces,[3 2 1]),numel(temp_ctx_trials),[]);

%                     Z = pdist(temp_ctx_traces_flat_trials, 'cosine');
%                     figure;
%                     imagesc(1-squareform(Z));
%                     title(sprintf('%s, file %d, Cell to cell similarity, %ss',cond_name,n_file, ops.context_name{n_cntxt}));
%                     xlabel('Cells');
%                     ylabel('Cells');
%                     figure;
%                     dendrogram(linkage(Z));
%                     title(sprintf('%s, file %d, Cell to cell similarity, %ss',cond_name,n_file, ops.context_name{n_cntxt}));
% 
%                     Z = pdist(temp_ctx_traces_flat_cells, 'cosine');
%                     figure;
%                     imagesc(1-squareform(Z));
%                     title(sprintf('%s, file %d, Trial to Trial similarity, %ss',cond_name,n_file, ops.context_name{n_cntxt}));
%                     xlabel('Trials');
%                     ylabel('Trials');
%                     figure;
%                     dendrogram(linkage(Z));
%                     title(sprintf('%s, file %d, Trial to Trial similarity, %ss',cond_name,n_file, ops.context_name{n_cntxt}));
% 
%                     
                    figure;
                    hold on;
                    for ii = 1:numel(intr_cells_index2)
                        p = plot(temp_ctx_traces_flat_trials(ii,:)+ii);
                    end
                    for ii = 1:numel(temp_ctx_trials)
                        temp_frame = ii*ops.win.no_base_window_size - ops.win.no_base_window_size + 1;
                        line([temp_frame temp_frame], [1 numel(intr_cells_index2)+1], 'color', [0.4 0.4 0.4]);
                    end
                    axis tight;
                    title([cond_name ' ' ops.context_name{n_cntxt} ' trials, dataset ' num2str(n_file)]);
                    ylabel('Responsive cells');
                    xlabel('trial sequence (frames)');
                    
                    
                    figure;
                    hold on;
                    for ii = 1:numel(intr_cells_index2)
                        p = plot(temp_ctx_traces_flat_trials(ii,:)+ii, 'k');
                        if context_resp_cells2(intr_cells_index2(ii))>0
                            p.Color = 'r';
                        end
                    end
                    for ii = 1:numel(temp_ctx_trials)
                        temp_frame = ii*ops.win.no_base_window_size - ops.win.no_base_window_size + 1;
                        line([temp_frame temp_frame], [1 numel(intr_cells_index2)+1], 'color', [0.4 0.4 0.4]);
                    end
                    axis tight;
                    title([cond_name ' ' ops.context_name{n_cntxt} ' trials, dataset ' num2str(n_file)]);
                    ylabel('Responsive cells');
                    xlabel('trial sequence (frames)');



    %                 SI_concat = similarity_index([temp_ctx_traces_data'],[temp_ctx_traces_data']);
    %                 figure;
    %                 imagesc(SI_concat);
    %                 title('Similarity Index');
    %                 xlabel(sprintf('%d V1 cells', size(temp_ctx_traces_data,1)));

    %                 SI_concat = similarity_index([temp_ctx_traces_flat_trials],[temp_ctx_traces_flat_trials]);
    %                 figure;
    %                 imagesc(SI_concat);
    %                 title(sprintf('Cell to cell Similarity Index, %ss', ops.context_name{n_cntxt}));
    %                 xlabel(sprintf('Cells'));
    %                 ylabel(sprintf('Cells'));
    %                 
    %                 corr_mat_trial = zeros(size(temp_ctx_traces,3));
    %                 for ii = 1:size(temp_ctx_traces,3)
    %                     for jj = 1:size(temp_ctx_traces,3)
    %                         corr_mat_trial(ii,jj) = sum(diag(squeeze(temp_ctx_traces(:,:,ii))*squeeze(temp_ctx_traces(:,:,jj))'));
    %                     end
    %                 end
    %                 corr_mat_trial = corr_mat_trial * 1./diag(corr_mat_trial)
    %                 
    %                 figure; imagesc(corr_mat_trial);
    %                 title(sprintf('Trial to trial similarity, %ss', ops.context_name{n_cntxt}));
    %                 xlabel('Trials');
    %                 ylabel('Trials');
    %                 
    %                 corr_mat_trial_dist = 1 - corr_mat_trial;
    %                 Z = linkage(corr_mat_trial_dist);
    %                 figure;
    %                 dendrogram(Z)
                end

                % now for combined contexts

                temp_ctx_trials = find(logical((temp_trial_type == temp_context_mmn(1,1)) + (temp_trial_type == temp_context_mmn(2,1)) + (temp_trial_type == temp_context_mmn(3,1))));
                temp_ctx_traces = temp_traces_all(intr_cells2,:,temp_ctx_trials);
                temp_ctx_traces2 = temp_traces_all(dataset_index == n_file,:,temp_ctx_trials);
                temp_ctx_traces3 = temp_ctx_traces2(sort_x,:,:);
                
                
                temp_ctx_traces_flat_trials = reshape(temp_ctx_traces3, numel(intr_cells_index2), []);
                temp_ctx_traces_flat_cells = reshape(permute(temp_ctx_traces3,[3 2 1]),numel(temp_ctx_trials),[]);

                Z = pdist(temp_ctx_traces_flat_trials, 'cosine');
                figure;
                imagesc(1-squareform(Z));
                title(sprintf('%s, file %d, Cell to cell similarity all ctx',cond_name,n_file));
                xlabel(sprintf('Cells, %d dev', sum(context_resp_cells2(dataset_index == n_file,3)==1)));
                ylabel('Cells');
                figure;
                dendrogram(linkage(Z,'average'), 1000,'ColorThreshold',0.65);
                title(sprintf('%s, file %d, Cell to cell similarity all ctx',cond_name,n_file));
                
                [~, ~, dend_order] = dendrogram(linkage(Z,'average'), 1000,'ColorThreshold',0.65);
                Z = pdist(temp_ctx_traces_flat_trials(dend_order,:), 'cosine');
                figure;
                imagesc(1-squareform(Z));
                axis image;
                title(sprintf('%s, file %d, sorted Cell to cell similarity all ctx',cond_name,n_file));
                
                
                Z = pdist(temp_ctx_traces_flat_cells, 'cosine');
                figure;
                imagesc(1-squareform(Z));
                title(sprintf('%s, file %d, Trial to Trial similarity all ctx',cond_name,n_file));
                xlabel(sprintf('Trials %dc, %dr %dd',dset_num_ctx_trials(n_file,1), dset_num_ctx_trials(n_file,2), dset_num_ctx_trials(n_file,3)));
                ylabel('Trials');
                figure;
                dendrogram(linkage(Z, 'average'), 1000,'ColorThreshold',0.55);
                title(sprintf('%s, file %d, Trial to Trial similarity all ctx',cond_name,n_file));
                
                [~,~, dend_order] = dendrogram(linkage(Z, 'average'), 1000,'ColorThreshold',0.55);
                Z = pdist(temp_ctx_traces_flat_cells(dend_order,:), 'cosine');
                figure;
                imagesc(1-squareform(Z));
                title(sprintf('%s, file %d, Trial to Trial similarity all ctx',cond_name,n_file));
                xlabel(sprintf('Trials %dc, %dr %dd',dset_num_ctx_trials(n_file,1), dset_num_ctx_trials(n_file,2), dset_num_ctx_trials(n_file,3)));
                ylabel('Trials');
            end
%             % and now for all traces
%             
%             intr_cells_index3 = find(logical(num_dataset == n_file));
%             temp_traces_dataset = temp_traces_all(intr_cells_index3,:,:);
% 
%             temp_traces_dataset_flat_trials = reshape(temp_traces_dataset, numel(intr_cells_index3), []);
%             temp_traces_dataset_flat_cells = reshape(permute(temp_traces_dataset,[3 2 1]),size(temp_traces_dataset,3),[]);
% 
%             Z = pdist(temp_traces_dataset_flat_trials, 'cosine');
%             figure;
%             imagesc(1-squareform(Z));
%             title(sprintf('Cell to cell similarity all trials'));
%             xlabel('Cells');
%             ylabel('Cells');
%             figure;
%             dendrogram(linkage(Z), 1000);
%             title(sprintf('Cell to cell similarity all trials'));
% 
%             Z = pdist(temp_traces_dataset_flat_cells, 'cosine');
%             figure;
%             imagesc(1-squareform(Z));
%             title(sprintf('Trial to Trial similarity all trials'));
%             xlabel('Trials');
%             ylabel('Trials');
%             dendrogram(linkage(Z), 1000);
%             title(sprintf('Trial to Trial similarity all trials'));
%             
% %           % how about cell to cell similarity on control trials
%             temp_cont_trials = find(logical(temp_trial_type(1:400)));
%             temp_traces_dataset_cont = temp_traces_all(intr_cells_index3,:,temp_cont_trials);
% 
%             temp_traces_dataset_flat_trials = reshape(temp_traces_dataset_cont, numel(intr_cells_index3), []);
%             temp_traces_dataset_flat_cells = reshape(permute(temp_traces_dataset_cont,[3 2 1]),size(temp_traces_dataset_cont,3),[]);
% 
%             Z = pdist(temp_traces_dataset_flat_trials, 'cosine');
%             figure;
%             imagesc(1-squareform(Z));
%             title(sprintf('Cell to cell similarity all controls'));
%             xlabel('Cells');
%             ylabel('Cells');
%             figure;
%             dendrogram(linkage(Z), 1000);
%             title(sprintf('Cell to cell similarity all controls'));
% 
%             Z = pdist(temp_traces_dataset_flat_cells, 'cosine');
%             figure;
%             imagesc(1-squareform(Z));
%             title(sprintf('Trial to Trial similarity all controls'));
%             xlabel('Trials');
%             ylabel('Trials');
%             dendrogram(linkage(Z), 1000);
%             title(sprintf('Trial to Trial similarity all controls'));


%             Z = pdist([1 2;1 3;5 6])
%             Z1 = squareform(Z)
%             
%             Z2 = linkage(Z)
%             
%             figure;
%             dendrogram(Z2)
%             
            
            
            %clear temp_ctx_trials temp_ctx_traces_data temp_frame;


            % context_mmn(n_cell,:,:) = [temp_MMN_freq(2), 10 + params.redundent_to_analyze, 19; ...
            %                           temp_MMN_freq(1), 19 + params.redundent_to_analyze, 28]';


        end

    end
end