function [figobj] = plot_community_gridlines(A,bilateral,community)
    % 
    % plot_community_gridlines(A, rois_bilateral, rois_community)
    % 
    % Given community labels and laterality labels 
    % 
    % By Manjari Narayan
    % 
    
    h =  imagesc(A);
    %h = heatmap(A,community,community);
    set(gcf,'Position',[299   132   779   628])
    palette = (brewermap(100,'PuOr')); % PRGn, RdYlBu, PuOr
    colormap(palette);
    Amax = max(A(:)); Amin = min(A(:));
    set(colorbar,'Limits',[-1*max(abs([Amax,Amin])) 1*max(abs([Amax,Amin]))])
    
    p = size(A,1);
    if(isempty(bilateral))
        bilateral = cat(2,ones(1,p/2),ones(1,p/2)*2)'; 
    end
    
    bilateral_no = grp2idx(bilateral);
    community_no = grp2idx(community);
    communities = unique(community_no);
    community_idx = []; start_idx = [];
    for bb=1:length(unique(bilateral_no))
        for cc=1:length(communities)
            start_idx = cat(1,start_idx,min(find(bilateral_no==bb & community_no==cc)));
            tmp_idx = ceil(median(find(bilateral_no==bb & community_no==cc)));
            community_idx = cat(1,community_idx,tmp_idx);
        end
    end
    
    % Plot Grid Lines
    for ii=2:length(start_idx)
        ll(ii) = line([0.5 p+0.5],[start_idx(ii)-0.5 start_idx(ii)-0.5],'LineWidth',1.5, 'Color', [0 0 0]);
        up(ii) = line([start_idx(ii)-0.5 start_idx(ii)-0.5],[0.5 p+0.5],'LineWidth',1.5, 'Color', [0 0 0]);
    end
    
    set(gca,'TickLength', [0.001 0.01])
    set(get(h,'Parent'), 'XTick', ...
    community_idx, 'XTickLabels', community(community_idx), ...
    'FontSize',15,'FontName','Fira Sans')
    set(get(h,'Parent'),'YTick', ...
    community_idx, 'YTickLabels', community(community_idx));
    set(get(h,'Parent'),'XTickLabelRotation',45, 'YTickLabelRotation',-30)
    hold on;
    [commX commY] = grid_communities(community_no .*(bilateral==1));
    comm(1) = line(commX,commY,'LineWidth',2.5,'Color',[0 0 0]);
    [bilX bilY] = grid_communities(community_no .*(bilateral==2));
    comm(2) = line(bilX,bilY,'LineWidth',2.5,'Color',[0 0 0]);

    % hlines = {};
    % hlines{1} =
    %
    figobj = {};
    figobj.h = h;
    figobj.ll = ll;
    figobj.up = up;
    figobj.comm = comm;

    
end

function [X Y] = grid_communities(c)
    
    nc = max(c);
    X = [];
    Y = [];
    for i = 1:nc
        ind = find(c == i);
        if ~isempty(ind)
            mn = min(ind) - 0.5;
            mx = max(ind) + 0.5;
            x = [mn mn mx mx mn NaN];
            y = [mn mx mx mn mn NaN];
            X = [X, x]; 
            Y = [Y, y];
        end
    end 
    
    
end