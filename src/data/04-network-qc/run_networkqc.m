function run_networkqc()
    
    opts = {};
    opts.atlas = 'Schaefer100_Yeo7';
    opts.atlas_size = 100;
    opts.file = fullfile(pwd,'src','data','04-network-qc','config','Schaefer100_Yeo7_labels.csv');
    opts.datadir = fullfile(getenv('PI_SCRATCH'), ...
                            'COMET', ...
                            'CausalConnectome', ...
                            'derivatives', ...
                            'fmriprep-fsl', ...
                            'denoiser', ...
                            opts.atlas);
    
    community = readtable(opts.file);
    opts.Ci = community.communityno;
	p = length(opts.Ci);
	opts.bilateral = [ones(p/2,1); ones(p/2,1)*2];
    
	tasks = ls(opts.datadir,'-1');
	datafiles = regexp(tasks,'\n','split');
	tasks = regexp(tasks,'.mat','split');
	tasks = regexprep(tasks,{'collect_roitimeseries_','\n','ses-d1_','ses-d2_'},{'','','',''});
	tasks = tasks(1:end-1);
	
	for tt=17:length(tasks)
		% save qc outputs in 
		opts.savedir = fullfile(pwd,'reports/qc/networks',tasks{tt});
		mkdir(opts.savedir)
		% Call plotting library
		addpath(genpath(fullfile(pwd,'..','export_fig')))
		addpath(genpath(fullfile(getenv('PI_HOME'),'software','matlab','matlab-library','lib','BrewerMap')))
	
	    data = load(fullfile(opts.datadir, ...
	                datafiles{tt}));
		savetablename = fullfile(opts.savedir,sprintf('cc_fmri_networkqc_task-%s',data.condition));
		metrics = {};
		metrics.subject = cell(length(data.X),1); 
		metrics.bilateral = nan(length(data.X),1); 
		metrics.ratio_cut = nan(length(data.X),1); 
		metrics.flow_cut = nan(length(data.X),1); 
	
		disp(['Starting ' tasks{tt} '...'])
	    for ii=1:length(data.X)
	        if(~isempty(data.X{ii}))
            
	            network = compute_correlation(data.X{ii},opts);
            
				% plot network matrix
				[figobj] = cc_fmri_plotmatrix(network,opts);
				savefilename = fullfile(opts.savedir,sprintf('cc_fmri_networkqc_task-%s_subject-%s',data.condition,data.subjects{ii}));
			
				% Collect bilateral and partitioning metrics
	            bilateral_metrics = compute_interhemispheric_agreement(network,opts);            
	            partitioning_metrics = compute_partitioning_score(network,opts);
            
	            disp(data.subjects{ii})
	            metrics_summary = sprintf('Bilateral: %2.3f, Ratio-cut: %2.3f, Conductance %2.3f',...
	                 bilateral_metrics.agreement, ....
	                 partitioning_metrics.ratio_cut, ....
	                 partitioning_metrics.conductance);
			
				title(metrics_summary); 
				disp(metrics_summary); 
				export_fig(savefilename,'-q98','-transparent');
            
				 metrics.subject{ii} = data.subjects{ii}; 
				 metrics.bilateral(ii) = round(bilateral_metrics.agreement,4);
				 metrics.ratio_cut(ii) = round(partitioning_metrics.ratio_cut,4);
				 metrics.flow_cut(ii) = round(partitioning_metrics.conductance,4);
				save(savefilename,'network','bilateral_metrics','partitioning_metrics');
				close all;
			else
				metrics.subject{ii} = data.subjects{ii};
	        end 
	    end
		metrics = struct2table(metrics);
		write(metrics,opts.savedir);
	end
end



function collect_timeseries()
    
    
    
    
end


function network = compute_correlation(X,opts)
    
    addpath('~/MATLAB/ggmClass')
    network = covariance.mle_sample_covariance(X);
    
end



function bilateral_metrics = compute_interhemispheric_agreement(network,opts)
    % Schaefer100_Yeo7 on GroupHCP_1200, agreement: 0.4382
    
    p = floor(size(network,1)/2);
    dsize = max(5,round(p/10));
    
    offdiagonal = ones(p,p)-tril(ones(p,p),-dsize)-triu(ones(p,p),dsize);
    interhem = network(1:p,p+1:2*p);
    
    bilateral_metrics = {};
    bilateral_metrics.agreement = ...
                 sum(sum(interhem.*offdiagonal))/sum(sum(offdiagonal));
    
end



function partition_metrics = compute_partitioning_score(network,opts)
    
    % Example on HCP1200, Schaefer100_Yeo7
    % >> ratio_cut: 5.6805
    % >> conductance: .638
    % >> partition_metrics.mean_module_matrix
    %
    % ans =
    %
    %   169.2519       NaN       NaN       NaN       NaN       NaN       NaN
    %    96.7860  123.7626       NaN       NaN       NaN       NaN       NaN
    %   104.1354   88.3283  131.1186       NaN       NaN       NaN       NaN
    %    71.3127   69.9570   77.7012   80.4060       NaN       NaN       NaN
    %    16.7712   13.6550   14.3854    9.2381    9.5874       NaN       NaN
    %    49.2534   34.2762   52.4467   44.2842   12.0262   83.7374       NaN
    %    94.3371   80.3473   66.9570   49.3680   32.7101   81.0018  247.1503
    
    
    partition_metrics = {};
    n_communities = length(unique(opts.Ci));
    
    cut_matrix = nan(n_communities,n_communities);
    size_matrix = nan(1,n_communities);
    ratio_matrix = cut_matrix;
    conductance_matrix = cut_matrix;
    for ii=1:n_communities
        size_matrix(ii) = sum(opts.Ci==ii);
        for jj=1:ii
            cut_matrix(ii,jj) = sum(sum(abs(network(opts.Ci==ii,opts.Ci==jj))));
            if(jj<ii)
                ratio_matrix(ii,jj) = cut_matrix(ii,jj)/(size_matrix(ii)*(size_matrix(jj)-1));
                conductance_matrix(ii,jj) = cut_matrix(ii,jj)/min(cut_matrix(ii,ii),cut_matrix(jj,jj));
            end
        end
    end
    %vol_matrix = diag(diag(cut_matrix));
    %cut_matrix = tril(cut_matrix,-1);
    
    partition_metrics.ratio_cut = nansum(nansum(ratio_matrix));
    partition_metrics.mean_module_matrix = cut_matrix;
    partition_metrics.conductance = ...
             nansum(nansum(conductance_matrix))/nchoosek(n_communities,2);

    
end
