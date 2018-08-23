function cc_fmri_subjects_conditions()
    % 
    % Returns subjects x conditions table
    % 
    DATADIR=['/Volumes/MACBOOKUSB/Datasets/' ...
                    'tms-fMRI/CC/roitimeseries'];
    %DATADIR = fullfile('data','interim','CC','roitimeseries');                

    
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
    end

    conditions.Properties.RowNames = regexprep(subjects_union',{'CausCon_'},'');
    
    writetable(conditions,'cc_fmri_subjects_conditions.csv','WriteVariableNames',1,'Delimiter',',');
    
end