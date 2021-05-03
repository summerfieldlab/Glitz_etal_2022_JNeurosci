%% This script is going to run an RSA at the time of outcome presentation
%  that is specified. It will split trials up into quartiles by their
%  probability of the presented door leading to the gain/loss state. Then
%  we will compare the representational similarity of identical and
%  non-identical transition probability bins across the two contexts in
%  each of our conditions (dependent and independent). Crossvalidation will
%  be performed by excluding the same session - same-session dissimilarity
%  quadrants when meaning across sessions.

%  c: Leonie Glitz, University of Oxford, 2020

% currently, the mask read in for this is the negative SPE mask, but you
% can change this manually

% TO-DO: adjust folder for data to read in once everything is on github

%% Step 1: load the masked voxel data

load(['/Users/leonieglitz/Desktop/Garrett&Glitz_etal_2021/outcome_fMRI_data/BOLD_outcome_SPEROI.mat']);
%this loads a struct
%BOLD(sub).sess(sess).gem(context).state(outcomeState).white_indices or .black_indices
%where data is a trial x voxel matrix with the outputs of GLM2a for the
%mask of interest, sorted by subject, session, context (==gem) and
%probability range (subjective probabilities of presented door transitioning
%to riksy state, divided into quartiles)

%% Step 2: compute the representational (dis)similarity between contexts and door chosen - outcome state combinations

emptyCounter = 0;

%loop through subjects
for sub = 1:29
    currentVoxelsReshaped =[];
 
    %loop through gems (contexts) and find corresponding beta.nii files
    for context = 1:4
        for outcomeState = 1:2
            for door = 1:2

                currentVoxelsLocal =[];

                %loop over sessions
                for sess = 1:4
                    currentVoxelsLocal = [];
                    BOLD_data = [];
                    
                    if door == 1
                        BOLD_data = BOLD(sub).sess(sess).gem(context).state(outcomeState).light_indices;
                    else
                        BOLD_data = BOLD(sub).sess(sess).gem(context).state(outcomeState).dark_indices;
                    end
                    
                    currentVoxelsLocal = BOLD_data;

                    if isempty(currentVoxelsLocal)
                        emptyCounter = emptyCounter +1;
                        currentVoxelsReshaped((context-1)*16+(outcomeState-1)*2 + door +(sess-1)*4,:) = NaN;
                    else
                        meanCurrentVoxels = nanmean(currentVoxelsLocal,1);
                        currentVoxelsReshaped((context-1)*16+(outcomeState-1)*2 + door +(sess-1)*4,:) = meanCurrentVoxels; %save as 16xvoxels matrix to be able to use pdist to obtain a 16x16 correlation matrix
                    end
                end
            end
        end
    end
    
    
    % going to use the correlation distance (1-correlation coefficient) to
    % compute dissimilarity
    %
    dissimVeryLarge(sub,:,:) = squareform(pdist(currentVoxelsReshaped,'Correlation'));
    
    dissimLargeDep(sub,:,:) = squareform(pdist(currentVoxelsReshaped(1:32,:),'Correlation'));
    dissimLargeDepBetween(sub,:,:) = dissimLargeDep(sub,1:16,17:32);
    simLargeDepSubBetween(sub,:,:) = 1- dissimLargeDepBetween(sub,:,:);
    
    dissimLargeIndep(sub,:,:) = squareform(pdist(currentVoxelsReshaped(33:64,:),'Correlation'));
    dissimLargeIndepBetween(sub,:,:) = dissimLargeIndep(sub,1:16,17:32);
    simLargeIndepSubBetween(sub,:,:) = 1-dissimLargeIndepBetween(sub,:,:);
    
end

['entries emtpy: ', num2str(emptyCounter)]

dissimLargeDep = dissimLargeDep([1:2,4:28],:,:);
dissimLargeDepBetween = dissimLargeDepBetween([1:2,4:28],:,:);
dissimLargeIndep = dissimLargeIndep([1:2,4:28],:,:);
dissimLargeIndepBetween = dissimLargeIndepBetween([1:2,4:28],:,:);
dissimVeryLarge = dissimVeryLarge([1:2,4:28],:,:);
simLargeDepSubBetween = simLargeDepSubBetween([1:2,4:28],:,:);
simLargeIndepSubBetween = simLargeIndepSubBetween([1:2,4:28],:,:);

%NaN out diagonal 4x4s (same block to same block) for crossvalidation and then mean
for i = 1:(64/4)
    dissimVeryLarge(:,1+((i-1)*4):(4+(i-1)*4),1+((i-1)*4):(4+(i-1)*4)) = NaN;
end

%% dependent condition dissimilarity extraction

%extract block subsections from NaNed out matrix - within is same context -
%same context (different sessions) and between is between contexts (diff
%sessions)

dissimLargeDepNaNDiag = dissimVeryLarge(:,1:32,1:32);
dissimLargeWithinNaNDiag = dissimLargeDepNaNDiag(:,1:16,1:16);
dissimLargeDepBetweenNaNDiag = dissimLargeDepNaNDiag(:,1:16,17:32);

%now we NaN out the diagonal 4x4s
dissimLargeDepBetweenNaNDiag(:,1:4, 1:4) = NaN;
dissimLargeDepBetweenNaNDiag(:,5:8, 5:8) = NaN;
dissimLargeDepBetweenNaNDiag(:,9:12, 9:12) = NaN;
dissimLargeDepBetweenNaNDiag(:,13:16, 13:16) = NaN; 
%% now the crossvalidation needs to be implemented for the dependent condition
%  to do this, we will first take the 4x4s along the diagonal (same run to
%  same run) and then average across the remaining 4x4s (crossvalidated, as
%  in the data for the similarities computed comes from different runs)


tmpWithinAcross = [];
tmpBetweenAcross =[];

for sub = 1:27
    for squaresDown = 1:4
        for squaresAcross = 1:4
            %within dependent block contexts
            tmpWithinAcross(sub, squaresAcross,squaresDown,:,:) = squeeze(dissimLargeWithinNaNDiag(sub,1+(squaresDown-1)*4:4*squaresDown,1+(squaresAcross-1)*4:4+(squaresAcross-1)*4));
            %between dependent block contexts
            tmpBetweenAcross(sub,squaresAcross,squaresDown,:,:) = squeeze(dissimLargeDepBetweenNaNDiag(sub,1+(squaresDown-1)*4:4*squaresDown,1+(squaresAcross-1)*4:4+(squaresAcross-1)*4));
            
        end
    end
    
end

for sub = 1:27
    %within dependent block contexts
    cvalMeanWithinDependentAcross(sub,:,:,:) = squeeze(nanmean(tmpWithinAcross(sub,:,:,:,:),2)); %mean across subs
    cvalMeanWithinDependent(sub,:,:) = squeeze(nanmean(cvalMeanWithinDependentAcross(sub,:,:,:),2));
    
    %across dependent block contexts
    cvalMeanBetweenDependentAcross(sub,:,:,:) = squeeze(nanmean(tmpBetweenAcross(sub,:,:,:,:),2)); %mean across subs
    cvalMeanBetweenDependent(sub,:,:) = squeeze(nanmean(cvalMeanBetweenDependentAcross(sub,:,:,:),2));
    
end

%% dependent - statistical testing

%transform dissimilarity to similarity (correlation)
cvalMeanWithinDepSimilarity = 1- cvalMeanWithinDependent;
cvalMeanBetweenDepSimilarity = 1 - cvalMeanBetweenDependent;

%Fisher transform values
cvalWithinDepSimFisher = atanh(cvalMeanWithinDepSimilarity);
cvalBetweenDepSimFisher = atanh(cvalMeanBetweenDepSimilarity);

%mean on and off-diagonal within context
cvalWithinDepMeanDiag = mean([cvalWithinDepSimFisher(:,1,1),cvalWithinDepSimFisher(:,2,2),cvalWithinDepSimFisher(:,3,3),cvalWithinDepSimFisher(:,4,4)],2);
cvalWithinDepMeanOffDiag = mean([cvalWithinDepSimFisher(:,1,2),cvalWithinDepSimFisher(:,1,3),cvalWithinDepSimFisher(:,1,4),cvalWithinDepSimFisher(:,2,1),...
    cvalWithinDepSimFisher(:,2,3),cvalWithinDepSimFisher(:,2,4),cvalWithinDepSimFisher(:,3,1),cvalWithinDepSimFisher(:,3,2),cvalWithinDepSimFisher(:,3,4),...
    cvalWithinDepSimFisher(:,4,1),cvalWithinDepSimFisher(:,4,2),cvalWithinDepSimFisher(:,4,3)],2);

%mean on and off-diagonal between contexts
cvalBetweenDepMeanDiag = mean([cvalBetweenDepSimFisher(:,1,1),cvalBetweenDepSimFisher(:,2,2),cvalBetweenDepSimFisher(:,3,3),cvalBetweenDepSimFisher(:,4,4)],2);
cvalBetweenDepMeanOffDiag = mean([cvalBetweenDepSimFisher(:,1,2),cvalBetweenDepSimFisher(:,1,3),cvalBetweenDepSimFisher(:,1,4),cvalBetweenDepSimFisher(:,2,1),...
    cvalBetweenDepSimFisher(:,2,3),cvalBetweenDepSimFisher(:,2,4),cvalBetweenDepSimFisher(:,3,1),cvalBetweenDepSimFisher(:,3,2),cvalBetweenDepSimFisher(:,3,4),...
    cvalBetweenDepSimFisher(:,4,1),cvalBetweenDepSimFisher(:,4,2),cvalBetweenDepSimFisher(:,4,3)],2);


%ttests  diagonal vs off-diagonal in the two conditions
[h,p,ci,tWithinGemDep] = ttest(cvalWithinDepMeanDiag,cvalWithinDepMeanOffDiag)
[h,p,ci,tBetweenGemDep] = ttest(cvalBetweenDepMeanDiag,cvalBetweenDepMeanOffDiag)

%calculate difference between diagonal and off-diagonal similarities between the two contexts in
%dependent condition (for later interaction analysis)
dependentBetweenDifference = cvalBetweenDepMeanDiag - cvalBetweenDepMeanOffDiag;

%% independent condition dissimilarity extraction

%extract block subsections from NaNed out matrix - within is same context -
%same context (different sessions) and between is between contexts (diff
%sessions)

dissimLargeIndepNaNDiag = dissimVeryLarge(:,33:64,33:64);
dissimLargeIndepWithinNaNDiag = dissimLargeIndepNaNDiag(:,1:16,1:16);
dissimLargeIndepBetweenNaNDiag = dissimLargeIndepNaNDiag(:,1:16,17:32);

%now we NaN out the diagonal 4x4s for the between context comparison
dissimLargeIndepBetweenNaNDiag(:,1:4, 1:4) = NaN;
dissimLargeIndepBetweenNaNDiag(:,5:8, 5:8) = NaN;
dissimLargeIndepBetweenNaNDiag(:,9:12, 9:12) = NaN;
dissimLargeIndepBetweenNaNDiag(:,13:16, 13:16) = NaN; 
%% now the crossvalidation needs to be implemented for the independent condition
%  to do this, we will first take the 4x4s along the diagonal (same run to
%  same run) and then average across the remaining 4x4s (crossvalidated, as
%  in the data for the similarities computed comes from different runs)

tmpWithinAcross = [];
tmpBetweenAcross =[];
for sub = 1:27
    for squaresDown = 1:4
        for squaresAcross = 1:4
            %within independent block contexts
            tmpWithinAcross(sub, squaresAcross,squaresDown,:,:) = squeeze(dissimLargeIndepWithinNaNDiag(sub,1+(squaresDown-1)*4:4*squaresDown,1+(squaresAcross-1)*4:4+(squaresAcross-1)*4));
            %between independent block contexts
            tmpBetweenAcross(sub,squaresAcross,squaresDown,:,:) = squeeze(dissimLargeIndepBetweenNaNDiag(sub,1+(squaresDown-1)*4:4*squaresDown,1+(squaresAcross-1)*4:4+(squaresAcross-1)*4));
            
        end
    end
    
end

for sub = 1:27
    %within independent block contexts
    cvalMeanWithinIndependentAcross(sub,:,:,:) = squeeze(nanmean(tmpWithinAcross(sub,:,:,:,:),2)); %mean across subs
    cvalMeanWithinIndependent(sub,:,:) = squeeze(nanmean(cvalMeanWithinIndependentAcross(sub,:,:,:),2));
    %across independent block contexts
    cvalMeanBetweenIndependentAcross(sub,:,:,:) = squeeze(nanmean(tmpBetweenAcross(sub,:,:,:,:),2)); %mean across subs
    cvalMeanBetweenIndependent(sub,:,:) = squeeze(nanmean(cvalMeanBetweenIndependentAcross(sub,:,:,:),2));
    
end

%% independent - statistical testing

%transform dissimilarity to similarity (correlation)
cvalMeanWithinIndepSimilarity = 1- cvalMeanWithinIndependent;
cvalMeanBetweenIndepSimilarity = 1 - cvalMeanBetweenIndependent;

%Fisher transform values
cvalWithinIndepSimFisher = atanh(cvalMeanWithinIndepSimilarity);
cvalBetweenIndepSimFisher = atanh(cvalMeanBetweenIndepSimilarity);

%mean on and off-diag within gem
cvalWithinIndepMeanDiag = mean([cvalWithinIndepSimFisher(:,1,1),cvalWithinIndepSimFisher(:,2,2),cvalWithinIndepSimFisher(:,3,3),cvalWithinIndepSimFisher(:,4,4)],2);
cvalWithinIndepMeanOffDiag = mean([cvalWithinIndepSimFisher(:,1,2),cvalWithinIndepSimFisher(:,1,3),cvalWithinIndepSimFisher(:,1,4),cvalWithinIndepSimFisher(:,2,1),...
    cvalWithinIndepSimFisher(:,2,3),cvalWithinIndepSimFisher(:,2,4),cvalWithinIndepSimFisher(:,3,1),cvalWithinIndepSimFisher(:,3,2),cvalWithinIndepSimFisher(:,3,4),...
    cvalWithinIndepSimFisher(:,4,1),cvalWithinIndepSimFisher(:,4,2),cvalWithinIndepSimFisher(:,4,3)],2);

%mean on and off-diag between gem
cvalBetweenIndepMeanDiag = mean([cvalBetweenIndepSimFisher(:,1,1),cvalBetweenIndepSimFisher(:,2,2),cvalBetweenIndepSimFisher(:,3,3),cvalBetweenIndepSimFisher(:,4,4)],2);
cvalBetweenIndepMeanOffDiag = mean([cvalBetweenIndepSimFisher(:,1,2),cvalBetweenIndepSimFisher(:,1,3),cvalBetweenIndepSimFisher(:,1,4),cvalBetweenIndepSimFisher(:,2,1),...
    cvalBetweenIndepSimFisher(:,2,3),cvalBetweenIndepSimFisher(:,2,4),cvalBetweenIndepSimFisher(:,3,1),cvalBetweenIndepSimFisher(:,3,2),cvalBetweenIndepSimFisher(:,3,4),...
    cvalBetweenIndepSimFisher(:,4,1),cvalBetweenIndepSimFisher(:,4,2),cvalBetweenIndepSimFisher(:,4,3)],2);


%ttests on similarity difference between diagonal and off-diagonal
[h,p,ci,tWithinGemIndep] = ttest(cvalWithinIndepMeanDiag,cvalWithinIndepMeanOffDiag)
[h,p,ci,tBetweenGemIndep] = ttest(cvalBetweenIndepMeanDiag,cvalBetweenIndepMeanOffDiag)

%calculate difference between diagonal and off-diagonal similarities between the two contexts in
%independent condition (for later interaction analysis)
independentBetweenDifference = cvalBetweenIndepMeanDiag - cvalBetweenIndepMeanOffDiag;

%% interaction analyses

%is there a difference in the difference between diagonal and off-diagonal
%between the two conditions
[h,p,ci,tDifferenceBtwBlocks] = ttest(dependentBetweenDifference,independentBetweenDifference)


%is there a difference in the absolute dissimilarity of the diagonals in
%the different contexts (in the two conditions)
[h,p,ci,tDifferenceBtwDiags] = ttest(cvalWithinDepMeanDiag,cvalWithinIndepMeanDiag)


%% are the same gems in crossval more similar than between-gems in crossval in dependent block?

%mean on and off-diag within gem
cvalWithinDepMeanDiag = mean([cvalWithinDepSimFisher(:,1,1),cvalWithinDepSimFisher(:,2,2),cvalWithinDepSimFisher(:,3,3),cvalWithinDepSimFisher(:,4,4)],2);
cvalWithinDepMeanOffDiag = mean([cvalWithinDepSimFisher(:,1,2),cvalWithinDepSimFisher(:,1,3),cvalWithinDepSimFisher(:,1,4),cvalWithinDepSimFisher(:,2,1),...
    cvalWithinDepSimFisher(:,2,3),cvalWithinDepSimFisher(:,2,4),cvalWithinDepSimFisher(:,3,1),cvalWithinDepSimFisher(:,3,2),cvalWithinDepSimFisher(:,3,4),...
    cvalWithinDepSimFisher(:,4,1),cvalWithinDepSimFisher(:,4,2),cvalWithinDepSimFisher(:,4,3)],2);

%mean on and off-diag between gem
cvalBetweenDepMeanDiag = mean([cvalBetweenDepSimFisher(:,1,1),cvalBetweenDepSimFisher(:,2,2),cvalBetweenDepSimFisher(:,3,3),cvalBetweenDepSimFisher(:,4,4)],2);
cvalBetweenDepMeanOffDiag = mean([cvalBetweenDepSimFisher(:,1,2),cvalBetweenDepSimFisher(:,1,3),cvalBetweenDepSimFisher(:,1,4),cvalBetweenDepSimFisher(:,2,1),...
    cvalBetweenDepSimFisher(:,2,3),cvalBetweenDepSimFisher(:,2,4),cvalBetweenDepSimFisher(:,3,1),cvalBetweenDepSimFisher(:,3,2),cvalBetweenDepSimFisher(:,3,4),...
    cvalBetweenDepSimFisher(:,4,1),cvalBetweenDepSimFisher(:,4,2),cvalBetweenDepSimFisher(:,4,3)],2);

[h,p,ci,differenceSameOtherDepStat] = ttest(cvalWithinDepMeanDiag-cvalWithinDepMeanOffDiag,cvalBetweenDepMeanDiag-cvalBetweenDepMeanOffDiag)