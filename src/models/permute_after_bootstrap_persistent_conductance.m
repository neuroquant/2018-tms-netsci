function [T Tperm results] = permute_after_bootstrap_persistent_conductance(C)
    % C is a clique/motif size x threshold x resamples x condition matrix
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
    
    community_tbl = readtable('Schaefer200_Yeo7_labels.csv');
    %communities = unique(community_tbl.communityno);
    %regexprep(communities,{'SomMot','DorsAttn','SalVentAttn','Cont','Default'},...
    %                        {'SMN','DAN','VAN','FPN','DMN'});
    communities = {'Vis','SMN','DAN','VAN','Limbic','FPN','DMN'};

    Cii=6;Cjj=7;
    metricname = [communities{Cii} '-' communities{Cjj}]
    metricfun = @(x)(x.metrics.conductances(:,:,Cii,Cjj)');
    
    methodname = 'corr';
    methodtype = fullfile(['networktype_' methodname],[methodname '_conductance']);
    [C labels] = collect_condition(methodtype,metricfun);
    results.C = C;
    labels2 = regexprep(labels,'_','.');   
    
    condNo1 = 5;
    condNo2 = 14;
    [maxT T] = conductance_test_statistic(C(:,:,:,condNo1),C(:,:,:,condNo2),1);
    maxT

    nperm = 2000;
    maxTperm = zeros(1,nperm);
    Tperm = zeros(size(T,1),size(T,2),nperm);
    tic;
    for permNo=1:nperm
        samples = permuter(2*size(C,3));
        nsamples = length(samples)/2;
        permC = cat(3,C(:,:,:,condNo1),C(:,:,:,condNo2));
        results.samples(permNo,:) = samples;
        [maxTperm(permNo) Tperm(:,:,permNo)] =  ...
                conductance_test_statistic(permC(:,:,samples(1:nsamples)), ...
                                           permC(:,:,samples(nsamples+1:2*nsamples)),1);
    end
    permtime = toc;
    sprintf('Testing took %6f seconds',permtime)
    %histogram(maxTperm,20); 
    assert(all(maxT==maxTperm)==0,'Permutations Not Working')
    
    pval = sum(abs(maxTperm)>abs(maxT))/nperm;
    sprintf('p-value for %s - %s:%.4f',labels2{condNo1},labels2{condNo2},pval)
    

    
    DATADIR=['/Volumes/MACBOOKUSB/Datasets/' ...
                    'tms-fMRI/CC/roitimeseries'];
    SAVEDIR = fullfile(DATADIR,'ggms',methodtype,metricname)
    mkdir(SAVEDIR)
    save(fullfile(SAVEDIR,['permstats_' labels{condNo1} '_' labels{condNo2}]),...
        'T','maxT','maxTperm','Tperm');
    
    figure;
    subplot(1,2,1); imagesc(T); caxis([0 1.1*maxT]); 
    colormap(brewermap(100,'PuOr'));
    colorbar('Location','SouthOutside');
    title(sprintf('p-value for %s - %s:%.4f',labels2{condNo1},labels2{condNo2},pval));
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

function [samples] = permuter(twosample)
    samples = randperm(twosample);
end


function [maxDiff Diff] = conductance_test_statistic(C1,C2,signDiff)
    
    C1 = permute(C1,[3 1 2]);
    C2 = permute(C2,[3 1 2]);
    if(isempty(signDiff)||signDiff>0)
        %disp('> Statistic')
        %Diff = (C1-C2).*(C1>C2);
        [~,~,~,stats] = ttest2(C1,C2,'tail','right','vartype','unequal');
        Diff = (squeeze(stats.tstat));
    elseif(signDiff<0)
        %disp('< Statistic')
        %Diff = (C2-C1).*(C2>C1);
        [~,~,~,stats] = ttest2(C1,C2,'tail','left','vartype','unequal');
        Diff = (squeeze(stats.tstat));
    elseif(signDiff==0)
        [~,~,~,stats] = ttest2(C1,C2,'tail','both','vartype','unequal');
        Diff = abs(squeeze(stats.tstat));
    end
    if(ndims(Diff)==3)
        Diff = mean(Diff,3)';
    end
    Diff(Diff==Inf) = NaN;
    maxDiff = nanmax(nanmax(Diff,[],1),[],2);
    
end


function C = get_all_conductances(resampled,metricfun)
    
    nresamples = length(resampled);
    C = [];
    condition_labels = {};
    
    for ii=1:nresamples
        C = cat(3,C,metricfun(resampled{ii}));
    end
    
end

function [C condition_labels] = collect_condition(methodtype,metricfun)
    % internal function: 
    % collect_condition(methodtype,metricfun)
    %
    % INPUT
    %
    % methodtype = fullfile(['networktype_' methodname],[methodname '_conductance'])
    
    DATADIR=['/Volumes/MACBOOKUSB/Datasets/' ...
                    'tms-fMRI/CC/roitimeseries'];
    LOADDIR=fullfile(DATADIR,'ggms',methodtype);

    tms_filenames = dir(fullfile(LOADDIR,'clique*.mat'));
    tms_filenames = {tms_filenames.name};

    condition_labels = {};
    C = [];
    for conditionNo=[3:5 13:15]%length(tms_filenames)
        tms_filename = fullfile(LOADDIR,tms_filenames{conditionNo});
        myfile = matfile(tms_filename);
        data.resampled = myfile.('resampled');
        if(~isfield(data,'resampled'))
            continue;
        end
        C(:,:,:,conditionNo) = get_all_conductances(data.resampled,metricfun);
        nresamples = length(data.resampled);
        condition_labels{conditionNo} = ...
             regexprep(tms_filenames{conditionNo},{'clique_conductance_','.mat'},{'',''});
        disp(condition_labels{conditionNo});    
        % condition_labels((conditionNo-1)*nresamples + 1:nresamples) = ...
        %      [repmat(regexprep({tms_filenames{conditionNo}},'clique_conductances_',''),[1 nresamples])];
    end
    %nanrows = mean(sum(isnan(C(:,:,:,1)),2),3)==size(C,2);
    %C = C(find(~nanrows),:,:,:);    
end