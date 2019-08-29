function [Shat tmsShat] = cc_fmri_tms_common_models()
    
    DATADIR=['/Volumes/MACBOOKUSB/Datasets/' ...
                    'tms-fMRI/CC/SchaeferYeo100/nuisance_4P+6P+6acompcor'];
    %DATADIR = fullfile('data','interim','CC','roitimeseries');    
    %DATADIR=fullfile('~/Downloads','partialcorr_conductance')

    methodname = 'corr';            
    SAVEDIR=fullfile(DATADIR,'ggms',['networktype_' methodname]);
    mkdir(SAVEDIR);
    ggm_results = struct();
    
    tms_filenames = dir(fullfile(SAVEDIR,'*ggms*.mat'));
    tms_filenames = {tms_filenames.name}
    nconditions = length(tms_filenames);

    tmsShat = [];
    tms_ni = [];
    ggm_results = struct();
    for conditionNo=1:nconditions
        tms_filename = fullfile(SAVEDIR,...
                        tms_filenames{conditionNo});
        load(tms_filename,'Shat'); 
        tmsShat(:,:,conditionNo) = mean(Shat,3);
        tms_ni(conditionNo) = size(Shat,3);
    end 
    tms_ni = reshape(tms_ni,[length(tms_ni) 1]);
    clear Shat;
    p = size(tmsShat,1);
    Shat = sum(tmsShat.*permute(repmat(tms_ni,[1 p p]),[2 3 1]),3);
    Shat = Shat/sum(tms_ni);
    
    description = {};
    description.tms_ni = 'subject sample size for each condition'; 
    description.Shat = 'MLE (appropriately weighted by sample size) across conditions';
    description.tmsShat = 'Group MLE Sample Correlation per Subject'
    
    % Save Results
    [~,tmpfilename] = fileparts(tms_filename);
    savefilename = 'common_tms_model';
    save(fullfile(SAVEDIR,savefilename),'-struct','ggm_results');   
    save(fullfile(SAVEDIR,savefilename),'tmsShat','Shat',...
                            'tms_ni','tms_filenames','description','-append'); 
    
    
end