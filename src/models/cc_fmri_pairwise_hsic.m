function [Sigma_ker Sigma_T Sigma_R] = cc_fmri_pairwise_hsic(X,subjectno)
    
    XTMS = standardize.successive_normalize(X{subjectno});
    
    other_subjects = setdiff([1:length(X)],subjectno)
    [Sigma_T Sigma_R] = helper_pairwise_hsic(X(other_subjects),10,5);
    
    Sigma = corr(XTMS);
    Sigma_ker = (XTMS'*Sigma_T*XTMS)/trace(Sigma_T);
    
end


function [Sigma_T Sigma_R] = helper_pairwise_hsic(X,bw_sigma1,bw_sigma2)
    % 
    % X is a cell array of non-empty data matrices
    % 
    
   [n p] = size(X{1}); 
   [m] = length(X);
   
   HSIC = zeros(n,n);
   HSIC2 = zeros(p,p);
   for ii=1:m
       ii
       XTMS1 = X{ii}; 
       for jj=ii+1:m
           XTMS2 = X{jj};
           [~,tmp_HSIC] = fasthsic(XTMS1,XTMS2);
           [~,tmp_HSIC2] = fasthsic(XTMS1',XTMS2',bw_sigma2,bw_sigma2);
           
           HSIC = HSIC + tmp_HSIC;
           HSIC2 = HSIC2 + tmp_HSIC2;
           
       end
   end
   HSIC = HSIC/nchoosek(m,2);
   HSIC2 = HSIC2/nchoosek(m,2);
   
   Sigma_T = HSIC;
   Sigma_R = HSIC2;
    
end