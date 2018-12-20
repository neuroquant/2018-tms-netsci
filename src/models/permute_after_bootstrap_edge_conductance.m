function [T Tperm results] = permute_after_bootstrap_edge_conductance(opts)
    % opts.Ci = [Cii Cjj]
    % opts.methodname 
    % opts.networktype
    % opts.metrictype 
    % opts.ATLAS 
    % (UNUSED) opts.C is a clique/motif size x threshold x resamples x condition matrix
    % 
    % - T is the observed test statistic for clique size x threshold x condition
    % - Tperm is the observed test statistic for clique size x threshold x nperms
    % 
    % Example: 
    % 
    s = RandStream('mt19937ar','Seed',1);
    RandStream.setGlobalStream(s);
    results = {};
    results.randstream = s;
    results.state = rng;

    if nargin==0
        Cii=6;Cjj=6;
        methodname = 'corr';
        networktype = 'partialcorr';
        metrictype = 'edge';
        ATLAS = 'SchaeferYeo100';
        
    elseif nargin>=1
        Cii=opts.Ci(1);
        Cjj=opts.Ci(2);
        methodname = opts.methodname;
        networktype = opts.networktype;
        metrictype = opts.metrictype;
        ATLAS = opts.ATLAS;
    end
    
    switch ATLAS
    case 'SchaeferYeo100'
        community_tbl = readtable('Schaefer100_Yeo7_labels.csv');
        %communities = unique(community_tbl.communityno);
        %regexprep(communities,{'SomMot','DorsAttn','SalVentAttn','Cont','Default'},...
        %                        {'SMN','DAN','VAN','FPN','DMN'});
    case 'SchaeferYeo200'
        community_tbl = readtable('Schaefer200_Yeo7_labels.csv');
    end
    communities = {'Vis','SMN','DAN','VAN','Limbic','FPN','DMN'};
    metricname = [communities{Cii} '-' communities{Cjj}]
    


    switch metrictype
    case 'clique'
        %% For persistent conductance
        metricfun = @(x)(x.metrics.conductances(:,:,Cii,Cjj)');
        methodtype = fullfile(['networktype_' methodname],[ networktype '_conductance'],'clique*');
        [C labels] = collect_condition(methodtype,metricfun,ATLAS);
        methodtype = fullfile(['networktype_' methodname],[ networktype '_conductance']);
        
    case 'edge'
        %% For edge conductance
        metricfun = @(x)(x.metrics(Cii,Cjj));
        methodtype = fullfile(['networktype_' methodname],[ networktype '_conductance'],...
        ['edge_conductance_*']);
        [C labels] = collect_condition(fullfile(methodtype),metricfun,ATLAS);
        methodtype = fullfile(['networktype_' methodname],[ networktype '_conductance'],...
        ['edge_conductance']);
        
    end

    DATADIR=fullfile(getenv('CC_DATADIR'),ATLAS,'roitimeseries');
    SAVEDIR = fullfile(DATADIR,'ggms',methodtype,metricname);
    mkdir(SAVEDIR)
    save([SAVEDIR filesep '/conductance_data.mat'],'C','labels','metrictype','methodtype','networktype');


    C = squeeze(C);
    results.C = C;
    labels = regexprep(labels,'edge_conductance_','');
    labels2 = regexprep(labels,'_','.');   
    nConditions = 16;
    
    tournamentG = zeros(length(labels),length(labels));
    for ii=1:nConditions
        for jj=setdiff(1:nConditions,ii)
            condNo1 = ii;
            condNo2 = jj;
            disp([labels{condNo1} '_' labels{condNo2}]);
            
            %%%%%%%%%%%%%%%%%%%%
            % DATADIR=['/Volumes/MACBOOKUSB/Datasets/' ...
            %                'tms-fMRI/CC/roitimeseries'];
            % LOADDIR = fullfile(DATADIR,'ggms',methodtype,metricname);
            %
            % switch metrictype
            % case 'edge'
            %     load(fullfile(LOADDIR,['edge_permstats_' labels{condNo1} '_' labels{condNo2}]));
            % case 'clique'
            %     load(fullfile(LOADDIR,['permstats_' labels{condNo1} '_' labels{condNo2}]));
            % end
            %%%%%%%%%%%%%%%%%%%
            if(ndims(C)==4)
            [maxT T] = ...
                 conductance_test_statistic(C(:,:,:,condNo1),C(:,:,:,condNo2),1);
            elseif(ndims(C)==2)
                 [maxT T] = ...
                      conductance_test_statistic(C(:,condNo1),C(:,condNo2),1);
            else
                error('Collection of metrics x condition failed');
            end

            nperm = 10000;
            maxTperm = zeros(1,nperm);
            Tperm = zeros(size(T,1),size(T,2),nperm);
            tic;
            for permNo=1:nperm
                if(ndims(C)==3)
                    samples = permuter(2*size(C,3));
                else
                    samples = permuter(2*size(C,1));
                end
                nsamples = length(samples)/2;
                results.samples(permNo,:) = samples;
                if(ndims(C)==4)
                    permC = cat(3,C(:,:,:,condNo1),C(:,:,:,condNo2));
                    [maxTperm(permNo) Tperm(:,:,permNo)] =  ...
                        conductance_test_statistic(permC(:,:,samples(1:nsamples)), ...
                                         permC(:,:,samples(nsamples+1:2*nsamples)),1);
                else
                    permC = cat(1,C(:,condNo1),C(:,condNo2));
                    [maxTperm(permNo) Tperm(permNo)] = ...
                         conductance_test_statistic(permC(samples(1:nsamples)), ...
                                         permC(samples(nsamples+1:2*nsamples)),1);
                end
            end
            permtime = toc;
            sprintf('Testing took %6f seconds',permtime)
            %histogram(maxTperm,20);
            assert(all(maxT==maxTperm)==0,'Permutations Not Working')

            pval = sum(abs(maxTperm)>abs(maxT))/nperm;
            sprintf('p-value for %s - %s:%.6f',labels2{condNo1},labels2{condNo2},pval)


            results.maxT = maxT;
            results.maxTperm = maxTperm;
            results.maxTpval = pval;
            results.T = T;
            if(ndims(Tperm)==3)
                results.Tperm = Tperm(:,:,1:1000);
                results.nullTperm = mean(Tperm,3);
                results.pval = sum(abs(Tperm)>repmat(abs(T),[1 1 size(Tperm,3)]),3)/nperm;
                results.pval(isnan(T)) = NaN;
            end
            switch metrictype
            case 'edge'
                save(fullfile(SAVEDIR,...
                ['edge_permstats_' labels{condNo1} '_' labels{condNo2}]),...
                        'T','maxT','maxTperm','pval');
            case 'clique'
                save(fullfile(SAVEDIR,...
                ['permstats_' labels{condNo1} '_' labels{condNo2}]),...
                        '-struct','results');
            end

            T = squeeze(T);
            Tperm = squeeze(Tperm);
            if(ndims(Tperm)==3)
                close all;
                figure;
                subplot(1,2,1); imagesc(T); caxis([0 1.1*maxT]);
                colormap(brewermap(100,'PuOr'));
                colorbar('Location','SouthOutside');
                title(sprintf('p-value for %s - %s:%.6f',labels2{condNo1},labels2{condNo2},pval));
                set(gca,'fontsize',12,'FontName','Fira Sans')
                subplot(1,2,2);
                imagesc(mean(Tperm,3)); caxis([0 1.1*maxT]);
                colormap(brewermap(100,'PuOr'));
                colorbar('Location','SouthOutside');
                %histogram(maxTperm,100,'Normalization','countdensity');
                title(sprintf('Permutation Null for %s - %s',labels2{condNo1},labels2{condNo2}));
                set(gca,'fontsize',12,'FontName','Fira Sans')
                savefig(fullfile(SAVEDIR,['permstats_' labels{condNo1} '_' labels{condNo2} '.fig']));
                export_fig(fullfile(SAVEDIR,['permstats_' labels{condNo1} '_' labels{condNo2} '.png']),...
                           '-transparent');
            end

            if(nanmin(pval(:))<.001)
                if(length(pval)>1)
                    tournamentG(ii,jj) = ...
                        max(max(T(pval>=0)));
                else
                    tournamentG(ii,jj) = maxT;
                end
                %tournamentG(jj,ii) = -maxT;
            else
                tournamentG(ii,jj) = 0;
            end
            
        end
    end 
    
    tournament_tbl = array2table(tournamentG);
    labels2
    tournament_tbl.Properties.VariableNames = labels;
    tournament_tbl.Properties.RowNames = labels;
    SAVEDIR = fullfile(DATADIR,'ggms',methodtype,metricname);
    writetable(tournament_tbl,...
        fullfile(SAVEDIR,'edge_graph_all_conditions.csv'),'WriteRowNames',1);
    
    
end

function [samples] = permuter(twosample)
    samples = randperm(twosample);
end


function [maxDiff Diff] = conductance_test_statistic(C1,C2,signDiff)
    
    if(ndims(C1)==3)
        C1 = permute(C1,[3 1 2]);
        C2 = permute(C2,[3 1 2]);
        stdC1 = squeeze(nanstd(C1,[],1));
        n1 = squeeze(sum(~isnan(C1),1));
        n2 = squeeze(sum(~isnan(C2),1));
        stdC2 = squeeze(nanstd(C2,[],1));
        nanmat = n1 <= 10 | n2 <= 10 | squeeze(stdC1<1e-10) | squeeze(stdC2<1e-10);
        stdC = (stdC1.*n1 + stdC2.*n2)./(n1 + n2 -2);
        stdC(nanmat==1) = 1.0;

        %% This is only a concern for the observed T-stat; not permuted null cases
        % if(~isempty(find(stdC<1e-12)))
        %     stdC(find(stdC<1e-12))
        %     warning('Numerical Errors in Conductance T-Statistic')
        % end
        
        % assert(isempty(find(stdC<1e-6)),'Numerical Errors in Conductance T-Statistic');
        % C1 = bsxfun(@rdivide,C1,stdC1);
        % C2 = bsxfun(@rdivide,C2,stdC2);
    else
        stdC1 = squeeze(nanstd(C1));
        n1 = squeeze(sum(~isnan(C1),1));
        n2 = squeeze(sum(~isnan(C2),1));
        stdC2 = squeeze(nanstd(C2));
        stdC = (stdC1.*n1 + stdC2.*n2)./(n1 + n2 -2);
    end


    if(isempty(signDiff)||signDiff>0)
        %disp('> Statistic')
        %Diff = squeeze((C1-C2).*(C1>C2));
        Diff = squeeze(nanmean((C1-C2).*(C1>C2),1))./stdC;
        if(ndims(C1)==3)
            Diff(nanmat==1) = NaN;
        end
        % C1(nanmat==1) = 0;
        % C2(nanmat==1) = 0;
        % [~,~,~,stats] = ttest2(C1,C2,'tail','right','vartype','unequal');
        % Diff = (squeeze(stats.tstat));
    elseif(signDiff<0)
        %disp('< Statistic')
        %Diff = (C2-C1).*(C2>C1);
        [~,~,~,stats] = ttest2(C1,C2,'tail','left','vartype','unequal');
        Diff = (squeeze(stats.tstat));
    elseif(signDiff==0)
        [~,~,~,stats] = ttest2(C1,C2,'tail','both','vartype','unequal');
        Diff = abs(squeeze(stats.tstat));
    end
    if(ndims(Diff)>=3)
        Diff = mean(Diff,3)';
    elseif(ndims(Diff)==2 & any(size(Diff)==1))
        Diff = mean(Diff,1);
    end
    Diff(Diff==Inf) = NaN;
    if(ndims(Diff)>=2)
        maxDiff = nanmax(nanmax(Diff,[],1),[],2);
    else
        maxDiff = Diff;
    end
end


function C = get_all_conductances(resampled,metricfun)
    
    nresamples = length(resampled);
    C = [];
    condition_labels = {};
    
    for ii=1:nresamples
        C = cat(3,C,metricfun(resampled{ii}));
    end
    
end

function [C condition_labels] = collect_condition(methodtype,metricfun,ATLAS)
    % internal function: 
    % collect_condition(methodtype,metricfun)
    %
    % INPUT
    %
    % methodtype = fullfile(['networktype_' methodname],[methodname '_conductance'])
    
    DATADIR=fullfile(getenv('CC_DATADIR'),ATLAS,'roitimeseries');

    LOADDIR=fullfile(DATADIR,'ggms',methodtype);

    tms_filenames = dir([LOADDIR '.mat']);
    tms_folder = fileparts(LOADDIR);
    tms_filenames = {tms_filenames.name};

    condition_labels = {};
    C = [];
    for conditionNo=1:length(tms_filenames)
        tms_filename = ...
             fullfile(tms_folder,tms_filenames{conditionNo});
        myfile = matfile(tms_filename);
        data.resampled = myfile.('resampled');
        if(~isfield(data,'resampled'))
            continue;
        end
        C(:,:,:,conditionNo) = get_all_conductances(data.resampled,metricfun);
        nresamples = length(data.resampled);
        condition_labels{conditionNo} = ...
             regexprep(tms_filenames{conditionNo},...
                 {'clique_conductance_','edge_conductance_','.mat'},...
                 {'','',''});
        disp(condition_labels{conditionNo});    
        % condition_labels((conditionNo-1)*nresamples + 1:nresamples) = ...
        %      [repmat(regexprep({tms_filenames{conditionNo}},'clique_conductances_',''),[1 nresamples])];
    end
    %nanrows = mean(sum(isnan(C(:,:,:,1)),2),3)==size(C,2);
    %C = C(find(~nanrows),:,:,:);    
end