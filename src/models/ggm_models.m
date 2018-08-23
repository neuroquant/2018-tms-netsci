function [ggm_results Shat] = ggm_models()
    
    DATADIR=['/Volumes/MACBOOKUSB/Datasets/' ...
                    'tms-fMRI/CC/roitimeseries'];
    %DATADIR = fullfile('data','interim','CC','roitimeseries');                
    SAVEDIR=fullfile(DATADIR,'ggms');
    mkdir(SAVEDIR);
    ggm_results = struct();
    
    tms_filenames = dir(fullfile(DATADIR,'*.mat'));
    tms_filenames = {tms_filenames.name};
    nconditions = length(tms_filenames);
    methodname = 'weightedcorr';
    SAVEDIR=fullfile(SAVEDIR,['networktype_' methodname]);
    mkdir(SAVEDIR);

    for conditionNo=1:nconditions
        tms_filename = fullfile(DATADIR,...
                        tms_filenames{conditionNo});
        Shat = get_cc_correlation_matrices(tms_filename,methodname); 
    
        ggm_results = struct();
        % tic;
        % ggm_results = estimator.model_average_populations(Shat);
        % ggm_time = toc;
        % disp(['Time Taken for Condition:' ggm_time])
        % ggm_results
    
        % Save Results
        [~,tmpfilename] = fileparts(tms_filename);
        savefilename = regexprep(tmpfilename,'collect_roitimeseries','stability_ggms');
        save(fullfile(SAVEDIR,savefilename),'-struct','ggm_results');
        save(fullfile(SAVEDIR,savefilename),'Shat','-append');
        
    end
    
    
end

function Shat = get_cc_correlation_matrices(filename,varargin)
    
    switch nargin
    case 2
        method = varargin{1};
    case 1
        method = 'corr';
    otherwise
        warning('More than 2 arguments not supported. Assuming first arg valid');
        method = 'corr';
    end
    
    studydata = load(filename)
    X = studydata.(studydata.Data);
    
    % Get subjects that are either NTHC or TEHC
    healthy_idx = find(cell2mat(cellfun(@(x)( ...
        ~isempty(strfind(x,'TEHC'))| ~isempty(strfind(x,'NTHC'))), ...
        studydata.subjects,...
        'UniformOutput',false))); 
    X = X(healthy_idx);
    [n p] = size(X{1});
    nsubjects = length(X);
    Shat = zeros(p,p,nsubjects);
    
    tms_covariate = readtable('cc_tms_ts_covariates.csv', ...
                    'ReadVariableNames',true,'Delimiter',',');
    Runc = {};
    if(isfield(tms_covariate,'onsets'))
        R = tms_covariate.ts(find(tms_covariate.onsets));
        Runc{1} = tms_covariate.onsets;
        canhrf = getcanonicalhrf(.4,2.4);
        Y = zeros(n,length(Runc));
        for ll=1:length(Runc)
            Y(:,ll) = conv([zeros(length(canhrf)-1,1); Runc{ll}],canhrf,'valid');
        end
    elseif(any(strcmp(tms_covariate.Properties.VariableNames,'ITI_1')))
        R = tms_covariate.ts(find(tms_covariate.ITI_1));
        Runc{1} = tms_covariate.ITI_1;
        Runc{2} = tms_covariate.ITI_2;
        Runc{3} = tms_covariate.ITI_3;
        Runc{4} = tms_covariate.ITI_4;
        Runc{5} = tms_covariate.ITI_5;
        Runc{6} = tms_covariate.ITI_6;
        canhrf = getcanonicalhrf(.4,2.4);
        Y = zeros(n,length(Runc));
        for ll=1:length(Runc)
            Y(:,ll) = conv([zeros(length(canhrf)-1,1); Runc{ll}],canhrf,'valid');
        end
    else
        warning('Could not find ITI information')
    end 
    
    
    for cc=1:nsubjects 
        Xnorm = standardize.successive_normalize(X{cc});
        switch method
        case 'kernelcorr'
            Shat(:,:,cc) = ...
                 kernel_correlation(X{cc},Xnorm,[Runc{:}]);
        case 'weightedcorr'
            Shat(:,:,cc) = weighted_correlation(Xnorm,Y);
        otherwise
            Shat(:,:,cc) = corr(Xnorm); 
        end
    end
end

function CorrXYX = weighted_correlation(Xnorm,Y)
    
    Yw = abs(Y); 
    Yw = bsxfun(@rdivide,Y,max(Y));
    [n p] = size(Xnorm);
    
    CorrXYX = zeros(p,p);
    for ll=1:size(Y,2)
        SigXYX = (Xnorm' * diag(Yw(:,ll)) * Xnorm);
        SigXYX = SigXYX/sum(Yw(:,ll));
        dSig = diag(SigXYX); 
        dSig(dSig<0) = 1;
        D_XYX = 1./sqrt(dSig);
        CorrXYX = CorrXYX + diag(D_XYX)*SigXYX*diag(D_XYX);
    end
    CorrXYX = CorrXYX/size(Y,2);
end
    

function KerCorrXYX = kernel_correlation(X,Xnorm,Y)
    
    % SigXY = sqSigXX * WXY * sqSigYY
    [KX,KY,KXbar,KYbar,~,HSIC] = fasthsic(Xnorm,Y);
    Wxy = (KXbar*KY*KXbar');
    KerSigXYX = Xnorm'*Wxy*Xnorm/(trace(Wxy));
    D_XYX = 1./sqrt(diag(KerSigXYX));
    KerCorrXYX = diag(D_XYX)*KerSigXYX*diag(D_XYX);
    
end


function add_paths()
    
    GGMDIR = fullfile(getenv('HOME'),'MATLAB','ggmClass');
    run(fullfile(GGMDIR,'setup.m'));
    
end