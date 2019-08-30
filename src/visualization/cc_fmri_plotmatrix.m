% Generic function to create interpretable connection matrices
function [figobj opts] = cc_fmri_plotmatrix(A,varargin)
	
    import utils.*
    
    opts = get_options();
    
    p = size(A,1); 
    
    figobj = plot_community_gridlines(A, ...
                                     opts.bilateral, ...
                                     opts.community.community);
end

function opts = get_options()
    
    opts = {};
    opts.atlas ='Schaefer100_Yeo7';
    opts.community = readtable(fullfile(pwd,'src','models','Schaefer100_Yeo7_labels.csv'));
    
    p = height(opts.community);
    opts.bilateral = [ones(p/2,1);ones(p/2,1)*2];
    
    
end