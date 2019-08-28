function cc_fmri_run_edge_conductance_tests()
    
    opts = {};
    opts.methodname = 'corr';
    opts.networktype = 'partialcorr';
    opts.metrictype = 'edge';
    opts.ATLAS = 'SchaeferYeo100';
    communities = {'Vis','SMN','DAN','VAN','Limbic','FPN','DMN'};
    ncommunities = length(communities);
    methodtype = fullfile(['networktype_' opts.methodname],[opts.networktype '_conductance'],...
        ['edge_conductance'])
    
    % for ii=1:ncommunities
    %     for jj=ii:ncommunities
    %         opts.Ci = [ii jj];
    %         permute_after_bootstrap_edge_conductance(opts);
    %     end
    % end
    
    collect_all_results = table();
    all_topsort_ranks = table();
    
    for ii=1:ncommunities
        for jj=ii:ncommunities
            topsort_ranks = [];
            metricname = [communities{ii} '-' communities{jj}]
            DATADIR=fullfile(getenv('CC_DATADIR'),opts.ATLAS,'roitimeseries');
            SAVEDIR = fullfile(DATADIR,'ggms',methodtype,metricname);
            results_filename = fullfile(SAVEDIR,'edge_graph_all_conditions.csv');
            tmp_table = readtable(results_filename,'ReadRowNames',1);
            tms_sites = tmp_table.Properties.VariableNames;
            G = digraph(table2array(tmp_table));
            toporder_idx = toposort(G,'order','stable');
            topsort_ranks(toporder_idx) = [length(toporder_idx):-1:1];
            all_topsort_ranks.(regexprep(metricname,'-','_')) = topsort_ranks';
            tmp_table.metricname = repmat({metricname},[height(tmp_table) 1]);
            tmp_table.Properties.RowNames = {};
            collect_all_results = vertcat(collect_all_results,tmp_table);
        end
    end
    
    all_topsort_ranks.Properties.RowNames = tms_sites;
    writetable(all_topsort_ranks,...
                fullfile('results','tournament_graphs', ...
                    opts.ATLAS,'roitimeseries','ggms',methodtype,'table_toporder_ranks.csv'),...
                'WriteRowNames',1);
    writetable(all_topsort_ranks,...
                fullfile('results','tournament_graphs', ...
                    opts.ATLAS,'roitimeseries','ggms',methodtype,'table_edge_graph_all.csv'),...
                'WriteRowNames',1);
    
end
