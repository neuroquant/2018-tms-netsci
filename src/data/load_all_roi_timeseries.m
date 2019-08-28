datadir = fullfile(getenv('PI_SCRATCH'),'COMET/CausalConnectome/derivatives/fmriprep-fsl/denoiser');
all_conditions = dir('00-bidsify/task-*.json');
all_conditions = regexprep({all_conditions.name},{'task-','_bold.json'},{'',''});

% % Schaefer100_Yeo7
% all_conditions = dir('00-bidsify/task-singlepulse*.json');
% all_conditions = regexprep({all_conditions.name},{'task-singlepulse','_bold.json'},{'',''});
%
% for conditionno=1:length(all_conditions)
%     condition = ['ses-d2_task-singlepulse' all_conditions{conditionno}]
%     load_bids_roi_timeseries(datadir,'condition',condition,'atlasname','Schaefer100_Yeo7');
% end
%
% condition = ['ses-d1_task-rest'];
% load_bids_roi_timeseries(datadir,'condition',condition,'atlasname','Schaefer100_Yeo7');

% Schaefer100_Yeo7

for conditionno=1:length(all_conditions)
    if(conditionno<=3)
        condition = ['ses-d1_task-' all_conditions{conditionno}];
    else
        condition = ['ses-d2_task-' all_conditions{conditionno}];
    end
    load_bids_roi_timeseries(datadir,'condition',condition,'atlasname','Schaefer100_Yeo7');
end

% Schaefer200_Yeo7
for conditionno=1:length(all_conditions)
    if(conditionno<=3)
        condition = ['ses-d1_task-' all_conditions{conditionno}];
    else
        condition = ['ses-d2_task-' all_conditions{conditionno}];
    end
    load_bids_roi_timeseries(datadir,'condition',condition,'atlasname','Schaefer200_Yeo7');
end

% Schaefer300_Yeo7
for conditionno=1:length(all_conditions)
    if(conditionno<=3)
        condition = ['ses-d1_task-' all_conditions{conditionno}];
    else
        condition = ['ses-d2_task-' all_conditions{conditionno}];
    end
    load_bids_roi_timeseries(datadir,'condition',condition,'atlasname','Schaefer300_Yeo7');
end

% Gordon333
for conditionno=1:length(all_conditions)
    if(conditionno<=3)
        condition = ['ses-d1_task-' all_conditions{conditionno}];
    else
        condition = ['ses-d2_task-' all_conditions{conditionno}];
    end
    load_bids_roi_timeseries(datadir,'condition',condition,'atlasname','Gordon333');
end

% Shen268
for conditionno=1:length(all_conditions)
    if(conditionno<=3)
        condition = ['ses-d1_task-' all_conditions{conditionno}];
    else
        condition = ['ses-d2_task-' all_conditions{conditionno}];
    end
    load_bids_roi_timeseries(datadir,'condition',condition,'atlasname','Shen268');
end

% Buckner7
for conditionno=1:length(all_conditions)
    if(conditionno<=3)
        condition = ['ses-d1_task-' all_conditions{conditionno}];
    else
        condition = ['ses-d2_task-' all_conditions{conditionno}];
    end
    load_bids_roi_timeseries(datadir,'condition',condition,'atlasname','Buckner7');
end

% Choi7
for conditionno=1:length(all_conditions)
    if(conditionno<=3)
        condition = ['ses-d1_task-' all_conditions{conditionno}];
    else
        condition = ['ses-d2_task-' all_conditions{conditionno}];
    end
    load_bids_roi_timeseries(datadir,'condition',condition,'atlasname','Choi7');
end