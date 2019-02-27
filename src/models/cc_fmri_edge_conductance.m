function [files] = cc_fmri_edge_conductance()
    % cc_fmri_edge_conductance
    % Computes only naive conductance at the edge level for all networks in DATADIR of the type specified by methodname
    % 
    % 
    % 
    ATLAS='Schaefer100_Yeo7';
    DATADIR=fullfile(getenv('CC_DATADIR'),ATLAS);
    
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
    switch ATLAS
    case {'Schaefer100_Yeo7', 'SchaeferYeo100'}
        community = readtable(['Schaefer100_Yeo7_labels.csv']);
    case {'Schaefer200_yeo7','SchaeferYeo200'}
        community = readtable(['Schaefer200_Yeo7_labels.csv']);
    end
    Ci = community.communityno;
    
    resampled = {};
    for conditionNo=1:length(tms_filenames)
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
             
             Asigned = A; (Asigned + Asigned')/2;         
             Apos = A.*(A>0); Apos = (Apos + Apos')/2; 
             Aneg = A.*(A<0); Aneg = (Aneg + Aneg')/2;
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
             
             Asigned = A; (Asigned + Asigned')/2;
             Apos = A.*(A>0); Apos = (Apos + Apos')/2; 
             Aneg = A.*(A<0); Aneg = (Aneg + Aneg')/2;
             Aabs = abs(A); A = (A + A')/2;
        end
        %thresholds = fliplr(unique(linspace(min_edge,max_edge,50)));
        
        [metrics] =  estimate_edge_conductance(A,Ci,0);
        [metrics_pos] =  estimate_edge_conductance(Asigned,Ci,1);
        [metrics_neg] =  estimate_edge_conductance(Asigned,Ci,-1);
        

         % Save Results
         [~,tmpfilename] = fileparts(tms_filename);
         savefilename = regexprep(tmpfilename,'stability_ggms','edge_conductance');
         save(fullfile(SAVEDIR,savefilename),'metrics','Apos', 'Aneg','metrics_pos','metrics_neg');
                
        for resampleNo=1:nresamples
            sprintf('Resample No: %d',resampleNo)
            if(use_partial_correlation)
                P = .5*(covariance.var_corr(graphs{resampleNo,1}(:,:,end)) + ...
                 covariance.var_corr(graphs{resampleNo,2}(:,:,end)));
                 Anew = (eye(size(P,1))-P);
                 Apos = Anew.*(Anew>0); Apos = (Apos + Apos')/2; 
                 Aneg = Anew.*(Anew<0); Aneg = (Aneg + Aneg')/2;
                 Asigned = Anew; (Asigned + Asigned')/2;
                 Anew = abs(Anew); Anew = (Anew + Anew')/2;  
             else
                 samples_idx = data.resampler.samples(resampleNo,:);
                 Anew = mean(data.Shat(:,:,samples_idx),3);
                 Apos = Anew.*(Anew>0); Apos = (Apos + Apos')/2; 
                 Aneg = Anew.*(Anew<0); Aneg = (Aneg + Aneg')/2;
                 Asigned = Anew; (Asigned + Asigned')/2;
                 Anew = abs(Anew); Anew = (Anew + Anew')/2;  
             end
            [resampled{resampleNo}.metrics] ...
                    = estimate_edge_conductance(Anew,Ci,0);
            [resampled{resampleNo}.metrics_pos] ...
                    = estimate_edge_conductance(Asigned,Ci,1);
            [resampled{resampleNo}.metrics_neg] ...
                    = estimate_edge_conductance(Asigned,Ci,-1);                
            if(mod(resampleNo,5)==0)
                save(fullfile(SAVEDIR,savefilename),'resampled','-append');
            end
        end
        save(fullfile(SAVEDIR,savefilename),'resampled','-append');

        clear resampled;
    end
    
end


function conductance = estimate_edge_conductance(A,Ci,signed)
% Calls community.cuts to compute conductance on the raw adjacency
    
    addpath('../../../netsci/netsci/')
    ncommunities = max(Ci); 
    conductance = zeros(ncommunities,ncommunities);
    
    for ii=1:ncommunities
        for jj=ii:ncommunities
            conductance(ii,jj) = ...
              community.cuts.conductance(A,Ci==ii,Ci==jj,[],signed);
        end
    end
    
    
end