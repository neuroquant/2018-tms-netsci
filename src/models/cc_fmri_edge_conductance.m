function [files] = cc_fmri_edge_conductance()
    % cc_fmri_edge_conductance
    % Computes only naive conductance at the edge level for all networks in DATADIR of the type specified by methodname
    % 
    % 
    % 
    
    DATADIR = getenv('DATADIR');
    
    methodname = 'corr';
    use_partial_correlation = true;
    if(use_partial_correlation)
        SAVEDIR=fullfile(DATADIR,'ggms',['networktype_' methodname],'partialcorr_conductance');
    else
        SAVEDIR=fullfile(DATADIR,'ggms',['networktype_' methodname],'corr_conductance');
    end
    mkdir(SAVEDIR);
    
    tms_filenames = dir(fullfile(DATADIR,'ggms',['networktype_' methodname],'*ggms*.mat'));
    tms_filenames = {tms_filenames.name};
    nconditions = length(tms_filenames);
    
    % Load Community Labels for ROIs
    community = readtable('Schaefer200_Yeo7_labels.csv');
    Ci = community.communityno;
    
    resampled = {};
    for conditionNo=1:15
        warning off
        tms_filename = fullfile(DATADIR, 'ggms',['networktype_' methodname], ...
                        tms_filenames{conditionNo});
        data = load(tms_filename);
    
        if(isfield(data,'stability_graphs'))
            Pi = data.stability_graphs(:,:,end);
        end
        if(use_partial_correlation)
            graphs = data.graphs;
            nresamples = size(graphs,1);
        else
            graphs = data.Shat;
            nresamples = size(data.resampler.samples,1);
        end
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
             edge_range = prctile(abs(edges),[50 98]);
             min_edge = round(edge_range(1),2);
             max_edge = round(edge_range(2),2);  
             A = abs(A); A = (A + A')/2;   
        else
             A = mean(data.Shat,3); edges = A(find(triu(A,1)));
             edge_range = prctile(abs(edges),[90 98]);
             min_edge = round(edge_range(1),2);
             max_edge = round(edge_range(2),2);
             if(sum(abs(A(:))>min_edge)/nchoosek(size(A,1),2)>.15)
                 disp('Updating Range')
                 edge_range = prctile(abs(edges),[95 99]);
                 min_edge = round(edge_range(1),2);
                 max_edge = round(edge_range(2),2);
             end
             A = abs(A); A = (A + A')/2;
        end
        %thresholds = fliplr(unique(linspace(min_edge,max_edge,50)));
        
        [metrics] = ...
             estimate_edge_conductance(A,Ci);;
             
         % Save Results
         [~,tmpfilename] = fileparts(tms_filename);
         savefilename = regexprep(tmpfilename,'stability_ggms','edge_conductance');
         save(fullfile(SAVEDIR,savefilename),'metrics');
                
        for resampleNo=1:nresamples
            sprintf('Resample No: %d',resampleNo)
            if(use_partial_correlation)
                P = .5*(covariance.var_corr(graphs{resampleNo,1}(:,:,end)) + ...
                 covariance.var_corr(graphs{resampleNo,2}(:,:,end)));
                 Anew = abs(eye(size(P,1))-P);
             else
                 samples_idx = data.resampler.samples(resampleNo,:);
                 Anew = mean(data.Shat(:,:,samples_idx),3);
                 Anew = abs(Anew);
             end
            [resampled{resampleNo}.metrics] ...
                    = estimate_edge_conductance(Anew,Ci);
            if(mod(resampleNo,5)==0)
                save(fullfile(SAVEDIR,savefilename),'resampled','-append');
            end
        end
        save(fullfile(SAVEDIR,savefilename),'resampled','-append');

        clear resampled;
    end
    
end


function conductance = estimate_edge_conductance(A,Ci)
% Calls community.cuts to compute conductance on the raw adjacency
    
    ncommunities = max(Ci); 
    conductance = zeros(ncommunities,ncommunities);
    
    for ii=1:ncommunities
        for jj=ii:ncommunities
            conductance(ii,jj) = community.cuts.conductance(A,Ci==ii,Ci==jj);
        end
    end
    
    
end