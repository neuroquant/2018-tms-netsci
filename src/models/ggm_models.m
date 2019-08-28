function [ggm_results Shat] = ggm_models()
    
    %DATADIR = fullfile('data','interim','CC','roitimeseries');
    DATADIR=fullfile(getenv('CC_DATADIR'),'SchaeferYeo100','roitimeseries');
    %DATADIR=fullfile(getenv('CC_DATADIR_ALT'),'Schaefer100_Yeo7');
                    
    SAVEDIR=fullfile(DATADIR,'ggms');
    mkdir(SAVEDIR);
    ggm_results = struct();
    
    tms_filenames = dir(fullfile(DATADIR,'*.mat'));
    tms_filenames = {tms_filenames.name};
    nconditions = length(tms_filenames);
    methodname = 'kernelcorr';
    SAVEDIR=fullfile(SAVEDIR,['networktype_' methodname]);
    mkdir(SAVEDIR);

    for conditionNo=1:nconditions
        tms_filename = fullfile(DATADIR,...
                        tms_filenames{conditionNo});
                        
        % Save Results
        [~,tmpfilename] = fileparts(tms_filename);
        savefilename = regexprep(tmpfilename,'collect_roitimeseries','stability_ggms');
        
        %ls(fullfile(SAVEDIR,[savefilename '.mat']))
        if(exist(fullfile(SAVEDIR,[savefilename '.mat'])))
            sprintf('Loading %s',fullfile(SAVEDIR,savefilename))
            load(fullfile(SAVEDIR,savefilename));
            
            A = mean(Shat,3); edges = A(find(triu(A,1)));
            edge_range = prctile(abs(edges),[80 98]);
            min_edge = round(edge_range(1),2);
            max_edge = round(edge_range(2),2);
            estimate_options = estimator.create_options();
            estimate_options.lambda_min = min_edge;
            estimate_options.lambda_max = max_edge;
            
        else
            Shat = get_cc_correlation_matrices(tms_filename,methodname);    
            ggm_results = struct();
            n_samples = size(Shat,3);
            estimate_options = estimator.create_options();
            resampling_options = estimate_options.resampler.options(n_samples);
            resampler = estimate_options.resampler.run(resampling_options);
            n_resamples = resampler.options.B;
            ggm_results.resampler = resampler;
            
            A = mean(Shat,3); edges = A(find(triu(A,1)));
            edge_range = prctile(abs(edges),[80 98]);
            min_edge = round(edge_range(1),2);
            max_edge = round(edge_range(2),2);
            estimate_options.lambda_min = min_edge;
            estimate_options.lambda_max = max_edge;
        end
        
        if(cond(mean(Shat,3))> 5e+3)
            for cc=1:size(Shat,3)
                alpha = .05;
                Shat(:,:,cc) = (1-alpha)*Shat(:,:,cc) + alpha*eye(size(Shat,1));
            end
            disp('Checking condition number')
            cond(mean(Shat,3))';
        end
       
        save(fullfile(SAVEDIR,savefilename),'Shat');
        
        
        % shrinkShat = zeros(size(Shat));
        % for cc=1:size(Shat,3)
        %     shrinkage = estimator.stein_shrinkage_sample_covariance(Shat(:,:,cc));
        %     tau_idx = max(find(shrinkage.tau>.001));
        %     shrinkShat(:,:,cc) = shrinkage.correlation{tau_idx};
        % end
        % shrinkage = estimator.stein_shrinkage_sample_covariance(mean(Shat,3))
        % tau_idx = max(find(shrinkage.tau>.001));
        % shrinkShat = shrinkage.correlation(tau_idx);
               
        estimate_options.path = estimate_options.lambdafun(min_edge,max_edge)
        tic;
        ggm_results = estimator.model_average_populations(Shat,estimate_options);
        ggm_time = toc;
        disp(['Time Taken for Condition:' ggm_time])
        ggm_results

        save(fullfile(SAVEDIR,savefilename),'-struct','ggm_results');
        save(fullfile(SAVEDIR,savefilename),'Shat','-append');
        
    end
    
    
end

function [Shat Xnorm] = get_cc_correlation_matrices(filename,varargin)
    
    switch nargin
    case 2
        method = varargin{1};
    case 1
        method = 'corr';
    otherwise
        warning('More than 2 arguments not supported. Assuming first arg valid');
        method = 'corr';
    end
    
    nTRs = 164;
    filename
    confounddir = fullfile(getenv('CC_DATADIR_ALT'),'..','confounds');
    studydata = load(filename);
    studydata.condition
    X = studydata.(studydata.Data);
    Xconfounds = studydata.confounds;
    
    % Get subjects that are either NTHC or TEHC
    notempty_idx = find(~cellfun(@isempty,X));
    healthy_idx = find(cell2mat(cellfun(@(x)( ...
        ~isempty(strfind(x,'TEHC'))| ~isempty(strfind(x,'NTHC'))), ...
        studydata.subjects,...
        'UniformOutput',false)));
    healthy_idx = intersect(healthy_idx,notempty_idx);
    X = X(healthy_idx);
    subjects = studydata.subjects(healthy_idx);
    if(~isempty(Xconfounds))
        Xconfounds = Xconfounds(healthy_idx);
    end
    
    [~, p] = size(X{1});
    nsubjects = length(X);
    Shat = zeros(p,p,nsubjects);
    
    tms_covariate = readtable('cc_tms_ts_covariates.csv', ...
                    'ReadVariableNames',true,'Delimiter',',');
    tms_covariate.onsets = ((tms_covariate.ITI>=1) & (tms_covariate.ITI<=2))*1.0;      
    tms_task = dlmread('cc_tms_task_covariate.csv');          
    % Runc = {};
    % if(isfield(tms_covariate,'onsets'))
    %     R = tms_covariate.ts(find(tms_covariate.onsets));
    %     Runc{1} = tms_covariate.onsets;
    %     canhrf = getcanonicalhrf(.4,2.4);
    %     Y = zeros(nTRs,length(Runc));
    %     for ll=1:length(Runc)
    %         Y(:,ll) = conv([zeros(length(canhrf)-1,1); Runc{ll}],canhrf,'valid');
    %     end
    % elseif(any(strcmp(tms_covariate.Properties.VariableNames,'ITI_1')))
    %     R = tms_covariate.ts(find(tms_covariate.ITI_1));
    %     Runc{1} = tms_covariate.ITI_1;
    %     Runc{2} = tms_covariate.ITI_2;
    %     Runc{3} = tms_covariate.ITI_3;
    %     Runc{4} = tms_covariate.ITI_4;
    %     Runc{5} = tms_covariate.ITI_5;
    %     Runc{6} = tms_covariate.ITI_6;
    %     canhrf = getcanonicalhrf(.4,2.4);
    %     Y = zeros(nTRs,length(Runc));
    %     for ll=1:length(Runc)
    %         Y(:,ll) = conv([zeros(length(canhrf)-1,1); Runc{ll}],canhrf,'valid');
    %     end
    % else
    %     warning('Could not find ITI information')
    % end
    
    
    for cc=1:nsubjects
        disp(sprintf('Subject %d',cc))
        subid = regexp(subjects{cc},'_','split');
        subid{2}
        confoundfile = dir(fullfile(confounddir, ...
                        ['*' subid{2} '*' ...
                        regexprep(studydata.condition,'_','') '*regressors.csv']));
        if(~isempty(confoundfile))
            confoundfile2 = fullfile(confounddir,confoundfile.name);
            fprep_confounds = readtable(confoundfile2,'Delimiter','\t', ...
                                        'TreatAsEmpty',{'n/a'});
            fprep_Xconfound = create_fprep_confounds(fprep_confounds, ...
                                                '4P+6M+6acompcor');    
        else
            fprep_Xconfound = [];
        end
        XX = X{cc}; 
        
        % Check that session length is exactly nTRs long or chose the last few TRs.
        if(size(XX,1)>nTRs)
            XX = XX(end-nTRs+1:end,:);
        end
        
        %% Used 11-2018
        % if(~isempty(Xconfounds))
        %     X_confound =  Xconfounds{cc}.X();
        %     X_perpY = Xy_orthogonalize(XX, X_confound(end-nTRs+1:end,4:9)); % Only 6 nuisance motion regressors
        %     Xnorm = standardize.successive_normalize(X_perpY);
        % else
        %     disp('Empty motion regressors')
        %     Xnorm = standardize.successive_normalize(XX);
        % end
        if(~isempty(fprep_Xconfound))
            fprep_Xconfound = zscore(fprep_Xconfound(end-nTRs+1:end,:));
            X_perpY = Xy_orthogonalize(XX,fprep_Xconfound); 
            Xnorm = standardize.successive_normalize(X_perpY);
        else
            disp('Empty motion regressors')
            Xnorm = standardize.successive_normalize(XX);
        end
        switch method
        case 'kernelcorr'
            Shat(:,:,cc) = ...
                 kernel_correlation(XX,Xnorm, cat(2,Y,Y.^2));
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


function X_perpY =  Xy_orthogonalize(X,Y)
% Nuisance signal regression
% X_perpY = Xy_orthogonalize(standardize.standardize_cols(X),options.nuisance);

	% p = size(X,2);
	% Xy = zeros(size(X));
	%
	% for ii=1:p
	% 	gsProj = gsr'*X(:,ii)/(gsr'*gsr)*gsr;
	% 	Xy = X(:,ii) - gsProj;
	% end
	
	mu = mean(X); 
	Xcen = bsxfun(@minus,X,mu); 
    
	proj = @(Z)(Z*pinv(Z'*Z)*Z');
	X_perpY = Xcen - proj(Y)*Xcen;
	X_perpY = bsxfun(@plus,X_perpY,mu);
 
end

function Xresiduals = nuisance_regression(X,Y,method)
    
    lmobj = fitlm(Y,X);
    Xresiduals = Y - X*lmobj.Coefficients;

end

function Xconfound = create_fprep_confounds(fprep_confounds,method)
    
    % Confounds List
    % {'csf'                       }
    % {'white_matter'              }
    % {'global_signal'             }
    % {'std_dvars'                 }
    % {'dvars'                     }
    % {'framewise_displacement'    }
    % {'t_comp_cor_00'             }
    % {'t_comp_cor_01'             }
    % {'t_comp_cor_02'             }
    % {'t_comp_cor_03'             }
    % {'t_comp_cor_04'             }
    % {'t_comp_cor_05'             }
    % {'a_comp_cor_00'             }
    % {'a_comp_cor_01'             }
    % {'a_comp_cor_02'             }
    % {'a_comp_cor_03'             }
    % {'a_comp_cor_04'             }
    % {'a_comp_cor_05'             }
    % {'cosine00'                  }
    % {'cosine01'                  }
    % {'cosine02'                  }
    % {'cosine03'                  }
    % {'cosine04'                  }
    % {'non_steady_state_outlier00'}
    % {'non_steady_state_outlier01'}
    % {'trans_x'                   }
    % {'trans_y'                   }
    % {'trans_z'                   }
    % {'rot_x'                     }
    % {'rot_y'                     }
    % {'rot_z'                     }
    % {'aroma_motion_01'           }
    % {'aroma_motion_02'           }
    % ....
    
    % Basic Confounds;
    % Default; '4P+6M'; % Discrete Cosines Basis for High Pass filtering 
                    % + 6 Movement Parameters + WM + CSF + dvars + FD

    
    Iscosines = strfind(fprep_confounds.Properties.VariableNames,'cosine');
    Iscosines = find(cellfun(@(x)(~isempty(x)),Iscosines,'UniformOutput',1));
    cosines = fprep_confounds.Properties.VariableNames(Iscosines);
    basic_confounds = fprep_confounds(:,...
                       [{'csf',...
                        'white_matter', ...
                        'std_dvars', ...
                        'framewise_displacement'} ...
                        cosines] ...                        
                        );
    motion =  fprep_confounds(:, { ...
                        'trans_x', ...
                        'trans_y', ...
                        'trans_z', ...
                        'rot_x', ...
                        'rot_y', ...
                        'rot_z' ...
                        });
    motion_squared = table();
    for ii=1:width(motion)
        colname = motion.Properties.VariableNames{ii};
        motion_squared.([colname '_sqd']) = motion.(colname).^2;
    end
    
    switch method
        
    case {'4P+6M+6acompcor','6acompcor'}
        additional_confounds = fprep_confounds(:,...
                       {'a_comp_cor_00',...
                        'a_comp_cor_01', ...
                        'a_comp_cor_02', ...
                        'a_comp_cor_03', ...
                        'a_comp_cor_04', ...
                        'a_comp_cor_05', ...
                        });
        additional_confounds = horzcat(additional_confounds, motion);                 
    case {'4P+12M+6acompcor'}
        additional_confounds = fprep_confounds(:,...
                       {'a_comp_cor_00',...
                        'a_comp_cor_01', ...
                        'a_comp_cor_02', ...
                        'a_comp_cor_03', ...
                        'a_comp_cor_04', ...
                        'a_comp_cor_05', ...
                        });
        additional_confounds = horzcat(additional_confounds, ...
                        motion, motion_squared);
                                        
    case {'4P+12M'}
        additional_confounds = horzcat(motion,motion_squared);        

    case {'4P+6M'}
        additional_confounds = motion;
        
    otherwise 
        disp('Unsupported argument')
        disp(method)
    end
    
    Xconfound = table2array(horzcat(basic_confounds, additional_confounds));
    
end


function add_paths()
    
    GGMDIR = fullfile(getenv('HOME'),'MATLAB','ggmClass');
    addpath(fullfile('GGMDIR','solvers','QUIC'));
    run(fullfile(GGMDIR,'setup.m'));
    
end