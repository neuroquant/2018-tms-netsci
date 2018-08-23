function resampled = cc_fmri_persistent_conductance()

    add_paths()
    
    DATADIR=['/Volumes/MACBOOKUSB/Datasets/' ...
                    'tms-fMRI/CC/roitimeseries/'];
    %DATADIR = fullfile('data','interim','CC','roitimeseries');                
    
    methodname = 'corr';
    use_partial_correlation = false;
    if(use_partial_correlation)
        SAVEDIR=fullfile(DATADIR,'ggms',['networktype_' methodname],'partialcorr_conductance');
    else
        SAVEDIR=fullfile(DATADIR,'ggms',['networktype_' methodname],'corr_conductance');
    end
    mkdir(SAVEDIR);
    
    tms_filenames = dir(fullfile(DATADIR,'ggms',['networktype_' methodname],'*.mat'));
    tms_filenames = {tms_filenames.name};
    nconditions = length(tms_filenames);
    
    % Load Community Labels for ROIs
    community = readtable('Schaefer200_Yeo7_labels.csv');
    Ci = community.communityno;
    
    
    for conditionNo=2:nconditions
        warning off
        tms_filename = fullfile(DATADIR, 'ggms',['networktype_' methodname], ...
                        tms_filenames{conditionNo});
        data = load(tms_filename);
    
        Pi = data.stability_graphs(:,:,end);
        graphs = data.graphs;
        % figure;
        % imagesc(A);
        % axis image off;
        % if(exist('brewermap'))
        %      colormap(flipud(brewermap([],'RdYlBu')));
        % end
        % colorbar;
        % title('Schaefer Yeo 100 Network on HCP1200 Group Data');
        % warning on;

        if(use_partial_correlation)
             P = covariance.var_corr(data.refit_estimator.inverse_covariance_estimate);
             A = eye(size(P,1))-P; edges = A(find(triu(A,1)));
             edge_range = prctile(abs(edges),[75 98]);
             min_edge = round(edge_range(1),2);
             max_edge = round(edge_range(2),2);  
             A = abs(A); A = (A + A')/2;   
        else
             A = mean(data.Shat,3); edges = A(find(triu(A,1)));
             edge_range = prctile(abs(edges),[90 98]);
             min_edge = round(edge_range(1),2);
             max_edge = round(edge_range(2),2);
             A = abs(A); A = (A + A')/2;
        end
        thresholds = fliplr(unique(linspace(min_edge,max_edge,50)));
        
        [metrics clique_adj cliques] = ...
             tda.persistent_conductance(A,Ci,thresholds,false,'persistent_conductance');
             
         % Save Results
         [~,tmpfilename] = fileparts(tms_filename);
         savefilename = regexprep(tmpfilename,'stability_ggms','clique_conductance');
         save(fullfile(SAVEDIR,savefilename),'metrics','clique_adj','cliques');
    
        % Plot DMN-Other Communities Conductance
        create_matrix_movie(squeeze(metrics.conductances(:,:,:,7)), ...
                            [fullfile(SAVEDIR,savefilename) '_DMN']); 
        create_matrix_movie(squeeze(metrics.conductances(:,:,1:6,6)), ...
                        [fullfile(SAVEDIR,savefilename) '_FPN']); 
                
        resampled = {};
        nresamples = size(graphs,1);

        for resampleNo=1:nresamples
            sprintf('Resample No: %d',resampleNo)
            if(use_partial_correlation)
                P = .5*(covariance.var_corr(graphs{resampleNo,1}(:,:,end)) + ...
                 covariance.var_corr(graphs{resampleNo,2}(:,:,end)));
                 Anew = abs(eye(size(P,1))-P);
             else
                 samples_idx = data.resampler.samples(resampleNo,:);
                 Anew = mean(data.Shat(:,:,samples_idx),3);
                 Anew = abs(A);
             end
            [resampled{resampleNo}.metrics,~,resampled{resampleNo}.cliques] ...
                = tda.persistent_conductance(Anew,Ci,thresholds,false,[]);
            if(mod(resampleNo,5)==0)
                save(fullfile(SAVEDIR,savefilename),'resampled','-append');
            end
        end
        save(fullfile(SAVEDIR,savefilename),'resampled','-append');
        
        clear resampled;
    end
    
end

function add_paths()
   
    NETSCIDIR = fullfile(getenv('HOME'),'MATLAB','netsci');
    addpath(NETSCIDIR) 
    addpath(fullfile(NETSCIDIR,'netsci'));
    addpath(fullfile(NETSCIDIR,'netsci','plot'));
    addpath(fullfile(NETSCIDIR,'netsci','utils'));
    addpath(fullfile(NETSCIDIR,'external','matlab-cliquer')); 
    run(fullfile(NETSCIDIR,'setup.m'));
    
    GGMDIR = fullfile(getenv('HOME'),'MATLAB','ggmClass');
    addpath(GGMDIR)
    
end