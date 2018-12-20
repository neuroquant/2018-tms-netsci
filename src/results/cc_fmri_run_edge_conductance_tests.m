function cc_fmri_run_edge_conductance_tests()
    
    env;
    opts = {};
    opts.methodname = 'corr';
    opts.networktype = 'partialcorr';
    opts.metrictype = 'edge';
    opts.ATLAS = 'SchaeferYeo100';
    
    communities = {'Vis','SMN','DAN','VAN','Limbic','FPN','DMN'};
    ncommunities = length(communities);
    
    for ii=1:ncommunities
        for jj=ii:ncommunities
            opts.Ci = [ii jj]; 
            permute_after_bootstrap_edge_conductance(opts);
        end
    end
end
