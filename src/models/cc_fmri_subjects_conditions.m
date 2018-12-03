function cc_fmri_subjects_conditions()
    % 
    % Returns subjects x conditions table
    % 
    DATADIR=fullfile(getenv('CC_DATADIR'),'SchaeferYeo100','roitimeseries');
    %DATADIR = fullfile('data','interim','CC','roitimeseries');                

    % Example
    % conditions = readtable('cc_fmri_subjects_conditions.csv');
    % imagesc(table2array(conditions))
    % colormap(viridis);
    % set(gca,'fontsize',16,'FontName','Fira Sans')
    % xlabel('Conditions'); ylabel('Subjects')
    %set(gca,'XTick',[1:14]+.5,'XTickLabels',regexprep(rowlabels,'_','.'),'XTickLabelRotation',45)

    tms_filenames = dir(fullfile(DATADIR,'*.mat'));
    tms_filenames = {tms_filenames.name};
    nconditions = length(tms_filenames);

    conditions = table();
     
    allsubjects = {};
    subjects_union = {};

    for conditionNo=1:nconditions
        tms_filename = fullfile(DATADIR,...
                        tms_filenames{conditionNo});
        load(tms_filename,'subjects')
        tms_subjects = regexp(subjects,'_tms_','split');
        for cc=1:length(tms_subjects)
            subjects{cc} = tms_subjects{cc}{1};
        end
        healthy_idx = find(cell2mat(cellfun(@(x)( ...
            ~isempty(strfind(x,'TEHC'))| ~isempty(strfind(x,'NTHC'))), ...
            subjects,...
            'UniformOutput',false))); 
        %length(healthy_idx);
        allsubjects{conditionNo} = subjects(healthy_idx);
        subjects_union = union(subjects_union,allsubjects{conditionNo});
        clear subjects
    end

    tms_labels = ...
         regexprep(tms_filenames,{'.mat','collect_roitimeseries_'},{'',''});
    for conditionNo=1:nconditions
        [isloa] = ismember(subjects_union,allsubjects{conditionNo});
        conditions.(tms_labels{conditionNo}) = isloa;
        clear isloa
    end

    conditions.Properties.RowNames = regexprep(subjects_union',{'CausCon_'},'');
    writetable(conditions,'cc_fmri_subjects_conditions.csv', ...
         'WriteVariableNames',1,'Delimiter',',','WriteRowNames',1);  
    
    figure;
    set(gcf,'Position',[440   148   776   650])     
    imagesc(table2array(conditions)); 
    axis image tight;
    set(gca,'XTick',1:width(conditions),...
            'XTickLabel',conditions.Properties.VariableNames,...
            'XTickLabelRotation',45,...
            'fontsize',16);
    
end