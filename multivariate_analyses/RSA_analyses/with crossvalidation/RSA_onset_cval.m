%% This script is going to run an RSA at the time of door onset  
%  that is specified. It will split trials up into quartiles by their 
%  probability of the presented door leading to the gain/loss state. Then 
%  we will compare the representational similarity of identical and 
%  non-identical transition probability bins across the two contexts in 
%  each of our conditions (dependent and independent). Crossvalidation will
%  be performed by excluding the same session - same-session dissimilarity
%  quadrants when meaning across sessions. 

%  c: Leonie Glitz, University of Oxford, 2020

% currently, the mask read in for this is the negative SPE mask, if you
% want to run this for any other brain region, please change the BOLD.mat
% file read in

% TO-DO: adjust folder for data to read in once everything is on github

%% Step 1: load the masked voxel data

load('/Users/leonieglitz/Desktop/Garrett&Glitz_etal_2021/onset_fMRI_data/BOLD_onset_SPEROI.mat');
%this loads a struct BOLD(sub).sess(sess).gem(context).prob(probRange).data
%where data is a trial x voxel matrix with the outputs of GLM2a for the
%mask of interest, sorted by subject, session, context (==gem) and
%probability range (subjective probabilities of presented door transitioning
%to riksy state, divided into quartiles)

%% Step 2: compute the representational (dis)similarity between different probability bins and contexts


subjectString = [1:20 22:27 29:31]; %this is important if you have exclusions


%loop over subjects
for sub = 1:length(subjectString)
    currentVoxelsReshaped =[];
    
    %loop over contexts and find corresponding beta.nii files
    for context = 1:4 
        %loop over probability ranges
        for probRange = 1:4
            currentBetas = [];
            currentVoxelsLocal =[];
            %loop over sessions
            for sess = 1:4
                currentVoxelsLocal = [];
                currentBetas = [];
                
                currentBetas = BOLD(sub).sess(sess).gem(context).prob(probRange).data;
                
                currentVoxelsLocal = currentBetas; 
 
                if isempty(currentVoxelsLocal)
                    'one empty condition';
                    currentVoxelsGlobal((probRange+(4*(sess-1)) + (context-1)*16),:) = repmat(NaN,([1,69]));
                else
                    currentVoxelsGlobal((probRange+(4*(sess-1)) + (context-1)*16),:) = nanmean(currentVoxelsLocal,1); %save sessions
                end
            end
        end 
    end
    
    
    % going to use the correlation distance (1-correlation coefficient) to
    % compute dissimilarity
   
    dissimVeryLarge(sub,:,:) = squareform(pdist(currentVoxelsGlobal,'Correlation'));
    
end

%NaN out diagonal 4x4s (same block to same block) 
for i = 1:(64/4)
    dissimVeryLarge(:,1+((i-1)*4):(4+(i-1)*4),1+((i-1)*4):(4+(i-1)*4)) = NaN;
end

%% dependent condition dissimilarity extraction

%extract block subsections from NaNed out matrix - within is same context -
%same context (different sessions) and between is between contexts (diff
%sessions)

dissimLargeDepNaNDiag = dissimVeryLarge(:,1:32,1:32); %all dependent block parts

dissimLargeDepWithinNaNDiag = dissimLargeDepNaNDiag(:,1:16,1:16); %same-to-same context dissimilarity Context 1
dissimLargeDepWithinNaNDiag2 = dissimLargeDepNaNDiag(:,17:32,17:32);%same-to-same context dissimilarity Context 2


dissimLargeDepBetweenNaNDiag = dissimLargeDepNaNDiag(:,1:16,17:32); %same-to-different context dissimilarity C1 - C2 (identical to C2-C1, which would be mirrored across diagonal)


%now we NaN out the diagonal 4x4s for the same
dissimLargeDepBetweenNaNDiag(:,1:4, 1:4) = NaN;
dissimLargeDepBetweenNaNDiag(:,5:8, 5:8) = NaN;
dissimLargeDepBetweenNaNDiag(:,9:12, 9:12) = NaN;
dissimLargeDepBetweenNaNDiag(:,13:16, 13:16) = NaN; 

%% now the crossvalidation needs to be implemented for the dependent condition
%  to do this, we will average across all 4x4s but the diagonal (crossvalidated, as
%  in the data for the similarities computed comes from different runs)

for sub = 1:29
    
    tmpDepWithin1 = squeeze(dissimLargeDepWithinNaNDiag(sub,:,:));
    tmpDepWithin2 = squeeze(dissimLargeDepWithinNaNDiag2(sub,:,:));
    
    %note: the same to same context RDMs are mirrored across the diagonal,
    %so we only need to include the lower diagonal 
    cvalMeanWithinDependent(sub,:,:) = nanmean(cat(3,...
        tmpDepWithin1(5:8, 1:4),...%context 1 to itself
        tmpDepWithin1(9:12, 1:4), tmpDepWithin1(9:12,5:8),...
        tmpDepWithin1(13:16,1:4), tmpDepWithin1(13:16,5:8), tmpDepWithin1(13:16,9:12),... 
        tmpDepWithin2(5:8, 1:4),... %now we move onto context 2
        tmpDepWithin2(9:12, 1:4), tmpDepWithin2(9:12,5:8),...
        tmpDepWithin2(13:16,1:4), tmpDepWithin2(13:16,5:8), tmpDepWithin2(13:16,9:12)),3);
    
    tmpDepBetween1 = squeeze(dissimLargeDepBetweenNaNDiag(sub,:,:));
    %because this is between contexts, the RDM is not mirrored across the
    %diagonal so we need to include the whole thing
    cvalMeanBetweenDependent(sub,:,:) = nanmean(cat(3,...
        tmpDepBetween1(1:4, 5:8), tmpDepBetween1(1:4, 9:12), tmpDepBetween1(1:4,13:16),... %context 1 to 2 (equivalent to 2-1, which would be mirrored)
        tmpDepBetween1(5:8, 1:4), tmpDepBetween1(5:8,9:12), tmpDepBetween1(5:8,13:16),...
        tmpDepBetween1(9:12, 1:4), tmpDepBetween1(9:12,5:8), tmpDepBetween1(9:12,13:16),...
        tmpDepBetween1(13:16,1:4), tmpDepBetween1(13:16,5:8), tmpDepBetween1(13:16,9:12)),3);

end




%% dependent - statistical testing

%transform dissimilarity to similarity (correlation)
cvalMeanWithinDepSimilarity = 1- cvalMeanWithinDependent;
cvalMeanBetweenDepSimilarity = 1 - cvalMeanBetweenDependent; 

%Fisher transform values
cvalWithinDepSimFisher = atanh(cvalMeanWithinDepSimilarity);
cvalBetweenDepSimFisher = atanh(cvalMeanBetweenDepSimilarity);

%replace any inf values with NaNs
cvalWithinDepSimFisher(isinf(cvalWithinDepSimFisher)) = NaN;
cvalBetweenDepSimFisher(isinf(cvalBetweenDepSimFisher)) = NaN;

%mean on and off-diag within context
cvalWithinDepMeanDiag = nanmean([cvalWithinDepSimFisher(:,1,1),cvalWithinDepSimFisher(:,2,2),cvalWithinDepSimFisher(:,3,3),cvalWithinDepSimFisher(:,4,4)],2);
cvalWithinDepMeanOffDiag = nanmean([cvalWithinDepSimFisher(:,1,2),cvalWithinDepSimFisher(:,1,3),cvalWithinDepSimFisher(:,1,4),cvalWithinDepSimFisher(:,2,1),...
    cvalWithinDepSimFisher(:,2,3),cvalWithinDepSimFisher(:,2,4),cvalWithinDepSimFisher(:,3,1),cvalWithinDepSimFisher(:,3,2),cvalWithinDepSimFisher(:,3,4),...
    cvalWithinDepSimFisher(:,4,1),cvalWithinDepSimFisher(:,4,2),cvalWithinDepSimFisher(:,4,3)],2);

%mean on and off-diag between contexts
cvalBetweenDepMeanDiag = nanmean([cvalBetweenDepSimFisher(:,1,1),cvalBetweenDepSimFisher(:,2,2),cvalBetweenDepSimFisher(:,3,3),cvalBetweenDepSimFisher(:,4,4)],2);
cvalBetweenDepMeanOffDiag = nanmean([cvalBetweenDepSimFisher(:,1,2),cvalBetweenDepSimFisher(:,1,3),cvalBetweenDepSimFisher(:,1,4),cvalBetweenDepSimFisher(:,2,1),...
    cvalBetweenDepSimFisher(:,2,3),cvalBetweenDepSimFisher(:,2,4),cvalBetweenDepSimFisher(:,3,1),cvalBetweenDepSimFisher(:,3,2),cvalBetweenDepSimFisher(:,3,4),...
    cvalBetweenDepSimFisher(:,4,1),cvalBetweenDepSimFisher(:,4,2),cvalBetweenDepSimFisher(:,4,3)],2);


%ttests
[h,p,ci,tWithinGemDep] = ttest(cvalWithinDepMeanDiag,cvalWithinDepMeanOffDiag)
[h,p,ci,tBetweenGemDep] = ttest(cvalBetweenDepMeanDiag,cvalBetweenDepMeanOffDiag)

%calculate difference between diagonal and off-diagonal similarities between the two contexts in
%dependent condition (for later interaction analysis)
dependentBetweenDifference = cvalBetweenDepMeanDiag - cvalBetweenDepMeanOffDiag;

%calculate difference between diagonal and off-diagonal similarities within
%a context (to check whether contexts are encoded consistently across
%conditions in crossvalidation)
dependentWithinDifference = cvalWithinDepMeanDiag - cvalWithinDepMeanOffDiag;

    
%% independentExtraction

%extract block subsections from NaNed out matrix
dissimLargeIndepNaNDiag = dissimVeryLarge(:,33:64,33:64);


dissimLargeIndepWithinNaNDiag = dissimLargeIndepNaNDiag(:,1:16,1:16); %same to same context 1
dissimLargeIndepWithinNaNDiag2 = dissimLargeIndepNaNDiag(:,17:32,17:32); %same to same context 2

dissimLargeIndepBetweenNaNDiag = dissimLargeIndepNaNDiag(:,1:16,17:32); % c1 - c2 (identical to c2-c1 just mirrored, so we do not need to include c2 -c1)

%now we NaN out the diagonal 4x4s for the between context comparison
dissimLargeIndepBetweenNaNDiag(:,1:4, 1:4) = NaN;
dissimLargeIndepBetweenNaNDiag(:,5:8, 5:8) = NaN;
dissimLargeIndepBetweenNaNDiag(:,9:12, 9:12) = NaN;
dissimLargeIndepBetweenNaNDiag(:,13:16, 13:16) = NaN; 


%% now the crossvalidation needs to be implemented for the independent condition
%  to do this, we will first take the 4x4s along the diagonal (same run to
%  same run) and then average across the remaining 4x4s (crossvalidated, as
%  in the data for the similarities computed comes from different runs)

for sub = 1:29
    
    tmpIndepWithin1 = squeeze(dissimLargeIndepWithinNaNDiag(sub,:,:));
    tmpIndepWithin2 = squeeze(dissimLargeIndepWithinNaNDiag2(sub,:,:));
    
    %note: the same to same context RDMs are mirrored across the diagonal,
    %so we only need to include the lower diagonal 
    
    cvalMeanWithinIndependent(sub,:,:) = nanmean(cat(3,...
        tmpIndepWithin1(5:8, 1:4),...  %context 1 to itself
        tmpIndepWithin1(9:12, 1:4), tmpIndepWithin1(9:12,5:8),...
        tmpIndepWithin1(13:16,1:4), tmpIndepWithin1(13:16,5:8), tmpIndepWithin1(13:16,9:12),... 
        ...
        tmpIndepWithin2(5:8, 1:4),... %now we move onto context 2
        tmpIndepWithin2(9:12, 1:4), tmpIndepWithin2(9:12,5:8), ...
        tmpIndepWithin2(13:16,1:4), tmpIndepWithin2(13:16,5:8), tmpIndepWithin2(13:16,9:12)),3);
    
    tmpIndepBetween1 = squeeze(dissimLargeIndepBetweenNaNDiag(sub,:,:));
    %because this is between contexts, the RDM is not mirrored across the
    %diagonal so we need to include the whole thing
    cvalMeanBetweenIndependent(sub,:,:) = nanmean(cat(3,...
        tmpIndepBetween1(1:4, 5:8), tmpIndepBetween1(1:4, 9:12), tmpIndepBetween1(1:4,13:16),... %context 1 to 2 (equivalent to 2-1, which would be mirrored)
        tmpIndepBetween1(5:8, 1:4), tmpIndepBetween1(5:8,9:12), tmpIndepBetween1(5:8,13:16),...
        tmpIndepBetween1(9:12, 1:4), tmpIndepBetween1(9:12,5:8), tmpIndepBetween1(9:12,13:16),...
        tmpIndepBetween1(13:16,1:4), tmpIndepBetween1(13:16,5:8), tmpIndepBetween1(13:16,9:12)),3);

end


%% independent - statistical testing

%transform dissimilarity to similarity (correlation)
cvalMeanWithinIndepSimilarity = 1- cvalMeanWithinIndependent;
cvalMeanBetweenIndepSimilarity = 1 - cvalMeanBetweenIndependent; 

%Fisher transform values
cvalWithinIndepSimFisher = atanh(cvalMeanWithinIndepSimilarity);
cvalBetweenIndepSimFisher = atanh(cvalMeanBetweenIndepSimilarity);


cvalWithinIndepSimFisher(isinf(cvalWithinIndepSimFisher)) = NaN;
cvalBetweenIndepSimFisher(isinf(cvalBetweenIndepSimFisher)) = NaN;

%mean on and off-diag within context
cvalWithinIndepMeanDiag = nanmean([cvalWithinIndepSimFisher(:,1,1),cvalWithinIndepSimFisher(:,2,2),cvalWithinIndepSimFisher(:,3,3),cvalWithinIndepSimFisher(:,4,4)],2);
cvalWithinIndepMeanOffDiag = nanmean([cvalWithinIndepSimFisher(:,1,2),cvalWithinIndepSimFisher(:,1,3),cvalWithinIndepSimFisher(:,1,4),cvalWithinIndepSimFisher(:,2,1),...
    cvalWithinIndepSimFisher(:,2,3),cvalWithinIndepSimFisher(:,2,4),cvalWithinIndepSimFisher(:,3,1),cvalWithinIndepSimFisher(:,3,2),cvalWithinIndepSimFisher(:,3,4),...
    cvalWithinIndepSimFisher(:,4,1),cvalWithinIndepSimFisher(:,4,2),cvalWithinIndepSimFisher(:,4,3)],2);

%mean on and off-diag between contexts
cvalBetweenIndepMeanDiag = nanmean([cvalBetweenIndepSimFisher(:,1,1),cvalBetweenIndepSimFisher(:,2,2),cvalBetweenIndepSimFisher(:,3,3),cvalBetweenIndepSimFisher(:,4,4)],2);
cvalBetweenIndepMeanOffDiag = nanmean([cvalBetweenIndepSimFisher(:,1,2),cvalBetweenIndepSimFisher(:,1,3),cvalBetweenIndepSimFisher(:,1,4),cvalBetweenIndepSimFisher(:,2,1),...
    cvalBetweenIndepSimFisher(:,2,3),cvalBetweenIndepSimFisher(:,2,4),cvalBetweenIndepSimFisher(:,3,1),cvalBetweenIndepSimFisher(:,3,2),cvalBetweenIndepSimFisher(:,3,4),...
    cvalBetweenIndepSimFisher(:,4,1),cvalBetweenIndepSimFisher(:,4,2),cvalBetweenIndepSimFisher(:,4,3)],2);


%ttests
[h,p,ci,tWithinGemIndep] = ttest(cvalWithinIndepMeanDiag,cvalWithinIndepMeanOffDiag)
[h,p,ci,tBetweenGemIndep] = ttest(cvalBetweenIndepMeanDiag,cvalBetweenIndepMeanOffDiag)

%calculate difference between diagonal and off-diagonal similarities between the two contexts in
%dependent condition (for later interaction analysis)
independentBetweenDifference = cvalBetweenIndepMeanDiag - cvalBetweenIndepMeanOffDiag;

%calculate difference between diagonal and off-diagonal similarities within
%a context (to check whether contexts are encoded consistently across
%conditions in crossvalidation)
independentWithinDifference = cvalWithinIndepMeanDiag - cvalWithinIndepMeanOffDiag;


%% interaction analyses

%is there a difference in the difference between diagonal and off-diagonal
%between the two conditions
[h,p,ci,tDifferenceBtwBlocks] = ttest(dependentBetweenDifference,independentBetweenDifference)

%is there a difference in the encoding of context identity between the two
%conditions (same to same) 
[h,p,ci,tDifferenceBtwBlocksWithin] = ttest(dependentWithinDifference,independentWithinDifference)


%is there a difference in the absolute similarity of the diagonals between conditions
[h,p,ci,tDifferenceBtwDiags] = ttest(cvalWithinDepMeanDiag,cvalWithinIndepMeanDiag)
