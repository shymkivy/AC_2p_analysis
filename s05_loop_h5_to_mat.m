clear;
close all;

addpath(genpath('C:\Users\ys2605\Desktop\stuff\caiman_sorter\caiman_sorter_functions'));

%%
%data_dir = 'F:\AC_data\caiman_data_dream3';
data_dir = 'F:\AC_data\caiman_data_missmatch';

ops = struct();
ops = f_cs_collect_ops_loop(ops);

%% evaluate params
ops.eval_params2.EvalSNRcaiman =            0;
ops.eval_params2.EvalSNR2 =                 1;
ops.eval_params2.EvalCNN =                  1;
ops.eval_params2.EvalRvalues =              0;
ops.eval_params2.EvalMinSigFrac =           0;
ops.eval_params2.EvalFiringStability =      0;

ops.eval_params2.RejThrSNRCaiman =          2;
ops.eval_params2.RejThrSNR2 =               5;
ops.eval_params2.RejThrCNN =                .97;
ops.eval_params2.RejThrRvalues =            0.5;
ops.eval_params2.RejThrMinSigFrac =         0.5;
ops.eval_params2.FiringStability =          0.01;

%% deconv params
ops.deconv.smooth_dfdt.params.convolve_gaus = 1;
ops.deconv.smooth_dfdt.params.gauss_kernel_simga = 100;
ops.deconv.smooth_dfdt.params.normalize = 0;
ops.deconv.smooth_dfdt.params.rectify = 1;

do_MCMC = 0;
ops.deconv.MCMC.params.save_SAMP = 0;

do_MCMC_cells = 'accepted'; % 'accepted', 'all'

%%
flist = dir([data_dir '\*.hdf5']);
fnames = {flist.name}';

for n_fl = 1:numel(fnames)
    [~, f_core, ext] = fileparts(fnames{n_fl});
    
    file_loc = [data_dir '\' fnames{n_fl}];
    
    f_sort_path = [data_dir '\' f_core '_sort.mat'];
    
    updates = 0;
    
    if exist(f_sort_path, 'file')
        fprintf('Updating %s... sort\n', f_core);
        load_data = load(f_sort_path);
        est = load_data.est;
        proc = load_data.proc;
        
        ops.eval_params_caiman = est.eval_params_caiman;
        ops.init_params_caiman = est.init_params_caiman;
    else
    
        fprintf('Creating new %s... sort\n', f_core);
        %% load data
        %if strcmpi(ext,'.hdf5') || strcmpi(ext,'.hdf') || strcmpi(ext,'.h5')
        temp_dims = h5read(file_loc,'/estimates/dims');

        est = f_cs_extract_h5_data(file_loc, temp_dims);
        est.dims = temp_dims;
    
        ops.eval_params_caiman = est.eval_params_caiman;
        ops.init_params_caiman = est.init_params_caiman;

        proc = f_cs_initialize_new_proc(est, ops);
        
        updates = 1;
    end 
    if ~updates
        flnames=fieldnames(ops.eval_params2);
        for n_field = 1:numel(flnames)
            if ops.eval_params2.(flnames{n_field}) ~= load_data.ops.eval_params2.(flnames{n_field})
                updates = 1;
            end
        end
        flnames=fieldnames(ops.deconv.smooth_dfdt.params);
        for n_field = 1:numel(flnames)
            if ops.deconv.smooth_dfdt.params.(flnames{n_field}) ~= load_data.ops.deconv.smooth_dfdt.params.(flnames{n_field})
                updates = 1;
            end
        end
    end
    
    %% evaluate comp
    if updates
        proc.comp_accepted_core = true(proc.num_cells,1);
        if ops.eval_params2.EvalSNRcaiman
            proc.comp_accepted_core = ((est.SNR_comp >= ops.eval_params2.RejThrSNRCaiman).*proc.comp_accepted_core);
        end
        if ops.eval_params2.EvalSNR2
            proc.comp_accepted_core = ((proc.SNR2_vals >= ops.eval_params2.RejThrSNR2).*proc.comp_accepted_core);
        end
        if ops.eval_params2.EvalCNN
            proc.comp_accepted_core = ((est.cnn_preds >= ops.eval_params2.RejThrCNN).*proc.comp_accepted_core);
        end
        if ops.eval_params2.EvalRvalues
            proc.comp_accepted_core = ((est.r_values >= ops.eval_params2.RejThrRvalues).*proc.comp_accepted_core);
        end
        if ops.eval_params2.EvalMinSigFrac
            proc.comp_accepted_core = (~(proc.num_zeros>=((1-ops.eval_params2.RejThrMinSigFrac)*proc.num_frames))).*proc.comp_accepted_core;
        end
        if ops.eval_params2.EvalFiringStability
            proc.comp_accepted_core = ((proc.firing_stab_vals >= ops.eval_params2.FiringStability).*proc.comp_accepted_core);
        end


        %% smoothdfdt
        proc.deconv.smooth_dfdt.S = zeros(proc.num_cells, proc.num_frames);
        smdfdt_params = ops.deconv.smooth_dfdt.params;

        data = double(est.C + est.YrA);
        fr = double(ops.init_params_caiman.data.fr);

        sigma_frames = smdfdt_params.gauss_kernel_simga*fr/1000;
        proc.deconv.smooth_dfdt.S = f_smooth_dfdt3(data, smdfdt_params.convolve_gaus, sigma_frames, smdfdt_params.normalize, smdfdt_params.rectify);

        for n_cell = 1:proc.num_cells
            temp_trace = proc.deconv.smooth_dfdt.S(n_cell,:);
            temp_trace(temp_trace <= 0) = [];
            proc.deconv.smooth_dfdt.S_std(n_cell) = sqrt(mean(temp_trace.^2));
        end
    end
    
    %% MCMC?
    if do_MCMC
        if ~isfield(proc.deconv, 'MCMC')
            proc.deconv.MCMC.S = cell(proc.num_cells,1);
            proc.deconv.MCMC.C = cell(proc.num_cells,1);
            proc.deconv.MCMC.SAMP = cell(proc.num_cells,1);
        else
            if ~isfield(proc.deconv.MCMC, 'S')
                proc.deconv.MCMC.S = cell(proc.num_cells,1);
            end
            if ~isfield(proc.deconv.MCMC, 'C')
                proc.deconv.MCMC.C = cell(proc.num_cells,1);
            end
            if ~isfield(proc.deconv.MCMC, 'SAMP')
                proc.deconv.MCMC.SAMP = cell(proc.num_cells,1);
            end
        end
        
        params.p = str2double(ops.deconv.MCMC.params.AR_val);
        params.f = double(ops.init_params_caiman.data.fr);
        params.Nsamples = ops.deconv.MCMC.params.Nsamples_param;
        params.B = ops.deconv.MCMC.params.B_param;

        dt = 1/double(ops.init_params_caiman.data.fr);

        f = waitbar(0,'Running MCMC...');
        for n_cell = 1:proc.num_cells

            if strcmpi(do_MCMC_cells, 'accepted')
                MCMC_list = proc.comp_accepted_core;
            elseif strcmpi(do_MCMC_cells, 'all')
                MCMC_list = true(proc.num_cells,1);
            end

            if MCMC_list(n_cell)
                if isempty(proc.deconv.MCMC.S{n_cell})
                    updates = 1;
                    params.sn = proc.noise(n_cell);

                    y = double(est.C(n_cell,:) + est.YrA(n_cell,:));

                    if params.p == 1
                        params.g = proc.gAR1(n_cell);
                    elseif params.p == 2
                        params.g = proc.gAR2(n_cell,:);
                    end

                    [SAMP, spikeRaster] = f_cs_compute_MCMC_core(y, params);

                    if ops.deconv.MCMC.params.save_SAMP
                        proc.deconv.MCMC.SAMP{n_cell} = SAMP;
                    end

                    if SAMP.process_ok
                        % need to fill in the zeros here for unprocessed signal
                        proc.deconv.MCMC.C{n_cell} = [zeros(1,SAMP.sig_start-1), mean(SAMP.C_rec,1)];
                        proc.deconv.MCMC.S{n_cell} = mean(spikeRaster,1);
                    end

                    waitbar(n_cell/proc.num_cells,f,'Running MCMC...');
                
                end
            end
        end
        close(f)
    end
    %%
    if updates
        disp('Saving...')
        save([data_dir '\' f_core '_sort.mat'], 'est', 'proc', 'ops', '-v7.3');
    else
        disp('No updates made, skipping...')
    end
end

disp('Done')
