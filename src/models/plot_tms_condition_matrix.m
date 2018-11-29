function plot_tms_condition_matrix(A,labels,filename)
    
    
    figure;
    set(gcf,'Position',[440   148   776   650])
    imagesc(A); 
    p = size(A,1);
    if(length(labels)~=p)
        error('Labels mismatched with matrix')
    end
    
    colorbar;
    colormap(brewermap(100,'PuOr'));
    
    set(gca,'XTick',1:p,'XTickLabel',labels,...
            'YTick',1:p,'YTickLabel',labels,...
            'XTickLabelRotation',45,...
            'YTickLabelRotation',45,...
            'fontsize',16,...
            'FontName','Fira Sans');
    
    export_fig([filename '.png'],'-q98','-transparent')