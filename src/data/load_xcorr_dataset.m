function obj = load_ts_dataset(filename,varargin)
%LOAD_TS_DATASET
% obj = load_ts_dataset(varargin)
% Loads a time-series dataset with expected variables outlined in a json configuration file.
% 
% USAGE: 
%   studydata = load_ts_dataset(filename,config)
%   
% INPUT
%   -filename
% 
%   (Optional)
%   - config contains fieldnames 'Data', 'Dimensions', 'PropertyNames'
% 
% DESCRIPTION: 
% Suppose you have stored a matrix of time-series in a variable called X 
% As long as template accurately describes the data variable as 
%   "Data": "X" 
% then the load_ts_dataset function will retrieve this information into the appropriate variable for subsequent analyses. 
% 
% Typically the data matrix has many dimensions such as timeseries, trials, subjects, channels or regions or voxels, etc..  In the absence of a uniform convention, the ordering of these variables should be specified in the variable `Dimensions`. 
% 
% The variable `PropertyNames` should contain the name of the expected string array, list or dictionary containing the names of the variables for each of the dimensions. For example, if the subject IDs of all subjects is stored in a variable `subID` and 3rd dimension of the data matrix corresponds to subjects then the third position of PropertyNames should contain "subID"
%   PropertyNames: "tslist", "roilist", "subID"
% 
% `Level1` and `Level2` explain the assumed nesting structure in the data for subsequent analyses. Loaded data will enforce the `Level1` measurements to take first dimension of the data matrix and `Level2` measurements to take the 3rd dimension. 
% 
% 2017, Manjari Narayan

    verbose = false;
    nargchk(1,2,nargin);
    if(nargin<2)
        cfg = jsonread('../../configs/sptms_eeg_xcorr_dataset_template.json');
    else
        try
            cfg = loadjson(varargin{1});
        catch me
            disp('Could not load json configuration file'); 
            disp(me);
        end            
    end
    templateVars = {'Data','Dimensions','PropertyNames'};
    
    fieldlist = fieldnames(cfg);
    DataVar = getfield(cfg,templateVars{1}); 
    Dimensions = getfield(cfg,templateVars{2}); 
    DataDims = length(Dimensions);
    PropertyVar = getfield(cfg,templateVars{3}); 
    assert(length(PropertyVar)==DataDims,'Length of PropertyNames does not match Dimensions');
      
    datastruct = load(filename);
    if(isfield(datastruct,'EEG'))
        datastruct = datastruct.EEG;
    end
    
    obj = struct();
    if(isfield(datastruct,DataVar))
        if(iscell(DataVar))
            for dd=1:length(DataVar)    
                obj = setfield(obj, ...
                        [templateVars{1} '_' num2str(dd)], ...
                        datastruct.(DataVar{dd}) ...
                        );
            end
        else
            obj = setfield(obj,templateVars{1},datastruct.(DataVar));
        end
    else
        obj = setfield(obj,templateVars{1},[]);
    end
    obj = setfield(obj,'Order',Dimensions);
    
    for dim_no=1:DataDims
        if(isfield(datastruct,PropertyVar{dim_no}))            
            obj = setfield(obj,PropertyVar{dim_no}, ...
                        datastruct.(PropertyVar{dim_no}));
        else
            obj = setfield(obj,PropertyVar{dim_no}, ...
                        []);
        end
    end 
    
    levelfields = sort(fieldlist( ...
                            cellfun(@(x)(~isempty(x)), ...
                                strfind(fieldlist,'Level')))); 
    
    function levelorder = get_levels(cfg,levelfields)
        levelorder = {};
        for levelno=1:length(levelfields)
            levelorder{levelno} = getfield(cfg,levelfields{levelno});
            if(verbose)
                disp([levelfields{levelno} ' is ' levelorder{levelno}]);
            end
            obj = setfield(obj,levelfields{levelno}, ...
                            getfield(cfg,levelfields{levelno}));
        end 
        levelorder = reshape(levelorder, [1 length(levelorder)]);       
    end    
    % for levelno=1:length(levelfields)
    %     obj = setfield(obj,levelfields{levelno}, ...
    %                     getfield(cfg,levelfields{levelno}));
    % end
    nolevels = any(cellfun(@isnan,get_levels(cfg,levelfields)));
    if(~nolevels)
        [locdim] = ismember(Dimensions, ...
                                get_levels(cfg,levelfields));
        assert(sum(locdim)==length(levelfields),'Mismatch between Dimensions and LevelX values')
    
        oldloc = 1:length(Dimensions);
        newloc = oldloc;
        if(verbose)
            disp(['Permute Data Matrix as ' num2str(newloc)]);  
        end
        obj = setfield(obj,templateVars{1}, ...
                               permute(obj.(templateVars{1}),newloc));
    end
end