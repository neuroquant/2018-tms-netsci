function results = load_bids_roi_timeseries(datadir,varargin)
    % Usage
    %   load_bids_roi_timeseries(datadir,'condition','ses-tp1_task-rest_run-2')
    % Inputs
    % datadir
    % subjectlist (optional): A cell of bids participant identifiers
    % atlasname : 'Schaefer200_Yeo7' (default), 'Schaefer100_Yeo7', 'Gordon333'
    % task : 'rest' (default) or another string.
    
    % $PI_SCRATCH/COMET/CausalConnectome/derivatives/denoiser/sub-xxx/sub-xxx_ses-d2_task-singlepulseRpMFG_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold_NR_roits
    p = inputParser;
    p.addRequired('datadir');
    p.addOptional('subjectlist',[]);
    p.addOptional('condition','ses-d2_task-singlepulseRpMFG');
    p.addOptional('atlasname','Schaefer100_Yeo7');
    
    p.parse(datadir,varargin{:});
    options = p.Results;
    subject_dir = dir([datadir filesep 'sub-*']);
    if(isempty(options.subjectlist))
        options.subjectlist = {subject_dir.name};
    else
        options.subjectlist = table2array(readtable(options.subjectlist,'ReadVariableName',0));
    end
    nsubjects = length(options.subjectlist);
    
    results = {};
    results.Data = 'X';
    results.X = {};
    results.subjects = {};
    results.files = {};
    results.confounds = {};
    results.condition = options.condition;
    taskname = regexp(options.condition,'_','split')
    results.task = regexprep(taskname(2),'task-','');
    
    for subjectNo=1:nsubjects
        
        %ses_dir = dir([datadir sprintf('/%s/*roits/',options.subjectlist{subjectNo})]);
        try
            results.subjects{subjectNo} = options.subjectlist{subjectNo};
            subjectfile = [datadir filesep options.subjectlist{subjectNo} ...
                            filesep '*' options.condition '*roits' ...
                            filesep '*' options.atlasname '*'];
            subjectfile = ls(subjectfile);
            subjectfile = subjectfile(1:end-1)
            % subjectfile = dir(subjectfile);
            % subjectfile = fullfile(subjectfile.folder,subjectfile.name);
            %disp(subjectfile)                
            results.files{subjectNo} = {subjectfile};
            tmpX = importdata(subjectfile);
            results.X{subjectNo} = tmpX;
        catch me
            results.files{subjectNo} = '';
            results.X{subjectNo} = [];
            disp(['Missing: ' options.subjectlist{subjectNo}]);
        end
    end
    disp(['Finished collecting timeseries for ' options.atlasname]);

    if(~exist([datadir filesep options.atlasname]))
        mkdir([datadir filesep options.atlasname]);
    end
    save([datadir filesep options.atlasname  ...
                  filesep 'collect_roitimeseries_' options.condition],'-struct','results');

end