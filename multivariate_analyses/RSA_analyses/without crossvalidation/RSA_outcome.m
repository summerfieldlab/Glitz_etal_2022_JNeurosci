function [pDep, statsDep, pIndep, statsIndep, pInteraction, statsInteraction] = RSA_outcome(maskData, fileName)
%% This function is going to run an RSA at the time of outcome in the ROI 
%  that is specified. It will split trials up by the door chosen (light/dark) 
%  and whether participants reached the gain/loss or the safe state. Then 
%  we will compare the representational similarity of identical and 
%  non-identical action-outcome pairings across the two contexts in 
%  each of our conditions (dependent and independent). 

%  cc: Leonie Glitz, University of Oxford, 2020

%  (a similar script for a predecessor of this analysis was written by Neil
%  Garrett too)

% TAKES INPUTS: maskData (i.e. mask path) and fileName -> save all data
% under what name 

%% Step 1: calculate beta indices of conditions of interest for each sub!

load('/Users/leonieglitz/Downloads/RW_flex_onset.mat');
behDat = readtable('/Volumes/Samsung_T5/gems/gem_dat.csv');

behDatNoNaN = behDat(~isnan(behDat.chooseLeft),:);
frequencyCount(1:29,1:16) = zeros; 

includedSubs = [1:20 22:27 29:31];

%we want to run separate RSAs for dependent and independent gems with
%themselves; within that we want to split up by door colour and outcome state


for sub = 1:29
    currSubject = behDatNoNaN(behDatNoNaN.participantID == includedSubs(sub),:);
    for sess = 1:4
        currSession = currSubject(currSubject.block_n == sess+1,:);
                
        indicesTmp = find(currSession.pick_black==0 & currSession.gem_presented==1 & currSession.outcomeState==1);
        subject(sub).sess(sess).gem(1).s(1).white_indices = indicesTmp;
        frequencyCount(sub,1) = frequencyCount(sub,1) + length(indicesTmp);
        
        indicesTmp = find(currSession.pick_black==1 & currSession.gem_presented==1 & currSession.outcome~=0);
        subject(sub).sess(sess).gem(1).s(1).black_indices = indicesTmp;
        frequencyCount(sub,2) = frequencyCount(sub,2) + length(indicesTmp);

        
        indicesTmp = find(currSession.pick_black~=1 & currSession.gem_presented==1 & currSession.outcome==0);
        subject(sub).sess(sess).gem(1).s(2).white_indices = indicesTmp;
        frequencyCount(sub,3) = frequencyCount(sub,3) + length(indicesTmp);
        
        indicesTmp = find(currSession.pick_black==1 & currSession.gem_presented==1 & currSession.outcome==0);
        subject(sub).sess(sess).gem(1).s(2).black_indices = indicesTmp;
        frequencyCount(sub,4) = frequencyCount(sub,4) + length(indicesTmp);
        
        
        
        indicesTmp = find(currSession.pick_black~=1 & currSession.gem_presented==2 & currSession.outcome~=0);
        subject(sub).sess(sess).gem(2).s(1).white_indices = indicesTmp;
        frequencyCount(sub,5) = frequencyCount(sub,5) + length(indicesTmp);

        
        indicesTmp = find(currSession.pick_black==1 & currSession.gem_presented==2 & currSession.outcome~=0);
        subject(sub).sess(sess).gem(2).s(1).black_indices = indicesTmp;
        frequencyCount(sub,6) = frequencyCount(sub,6) + length(indicesTmp);
        
        
        indicesTmp = find(currSession.pick_black~=1 & currSession.gem_presented==2 & currSession.outcome==0);
        subject(sub).sess(sess).gem(2).s(2).white_indices = indicesTmp;
        frequencyCount(sub,7) = frequencyCount(sub,7) + length(indicesTmp);
        
        
        indicesTmp = find(currSession.pick_black==1 & currSession.gem_presented==2 & currSession.outcome==0);
        subject(sub).sess(sess).gem(2).s(2).black_indices = indicesTmp;
        frequencyCount(sub,8) = frequencyCount(sub,8) + length(indicesTmp);
        
        
        
        
        indicesTmp = find(currSession.pick_black~=1 & currSession.gem_presented==3 & currSession.outcome~=0);
        subject(sub).sess(sess).gem(3).s(1).white_indices = indicesTmp;
        frequencyCount(sub,9) = frequencyCount(sub,9) + length(indicesTmp);

        
        indicesTmp = find(currSession.pick_black==1  & currSession.gem_presented==3 & currSession.outcome~=0);
        subject(sub).sess(sess).gem(3).s(1).black_indices = indicesTmp;
        frequencyCount(sub,10) = frequencyCount(sub,10) + length(indicesTmp);
        
        indicesTmp = find(currSession.pick_black~=1  & currSession.gem_presented==3 & currSession.outcome==0);
        subject(sub).sess(sess).gem(3).s(2).white_indices = indicesTmp;
        frequencyCount(sub,11) = frequencyCount(sub,11) + length(indicesTmp);
        
        indicesTmp = find(currSession.pick_black==1 & currSession.gem_presented==3 & currSession.outcome==0);
        subject(sub).sess(sess).gem(3).s(2).black_indices = indicesTmp;
        frequencyCount(sub,12) = frequencyCount(sub,12) + length(indicesTmp);
        
        
        indicesTmp = find(currSession.pick_black~=1 & currSession.gem_presented==4 & currSession.outcome~=0)
        subject(sub).sess(sess).gem(4).s(1).white_indices = indicesTmp;
        frequencyCount(sub,13) = frequencyCount(sub,13) + length(indicesTmp);
       
        indicesTmp = find(currSession.pick_black==1 & currSession.gem_presented==4 & currSession.outcome~=0);
        subject(sub).sess(sess).gem(4).s(1).black_indices = indicesTmp;
        frequencyCount(sub,14) = frequencyCount(sub,14) + length(indicesTmp);
        
        indicesTmp = find(currSession.pick_black~=1 & currSession.gem_presented==4 & currSession.outcome==0);
        subject(sub).sess(sess).gem(4).s(2).white_indices = indicesTmp;
        frequencyCount(sub,15) = frequencyCount(sub,15) + length(indicesTmp);
        
        indicesTmp = find(currSession.pick_black==1 & currSession.gem_presented==4 & currSession.outcome==0);
        subject(sub).sess(sess).gem(4).s(2).black_indices = indicesTmp;
        frequencyCount(sub,16) = frequencyCount(sub,16) + length(indicesTmp);
        
    end
    
 
    
end


%% Step 2: read in the appropriate beta files and compute the representational (dis)similarity


fs = filesep; %so script can be used on mac or windows

%setup
currentModelFolder = '/Volumes/Samsung_T5/gems/singleTrialModelNeil/Outcome'; %where your models are generally

subjectString = [1:20 22:27 29:31]; %this is important if you have exclusions

fileBase = 'beta_0'; %base for reading in the beta images
fileEnding = '.nii'; %ending
subPrefix = 'sub';

currentVoxelsReshaped = [];

relevantVoxels = struct(); %ensures that you don't get a struct error later
relevantVoxels.mask = spm_read_vols(spm_vol(maskData)); %read in mask for ROI - make sure same dimensions/resolution as your fMRI images
relevantVoxels.ROI = find(relevantVoxels.mask>0); %find where ROI is in the mask - these are the betas from the  ROI we want



%loop through subjects
for sub = 1:length(subjectString)
    currentVoxelsReshaped =[];
    currentVoxelsGlobal =[];
    %load the SPM.mat file
    load([currentModelFolder,...
        fs,subPrefix,num2str(subjectString(sub)),fs,'SPM.mat']); %load in the current subject's SPM file
    
    %loop through conditions/gems and find corresponding beta.nii files
        for gem = 1:4 %conditions
            for outcomeState = 1:2
                for door = 1:2
                    
                    currentBetas = [];
                    currentVoxels = [];
                    currentVoxelsLocal =[];
                    currentVoxelsGlobal =[]; 
                    for sess = 1:4
                        uniqueName = ['Sn(',num2str(sess),') outcome_onset_'];
                        
                        if door == 1
                            conditionCurr = subject(sub).sess(sess).gem(gem).s(outcomeState).white_indices;
                        else
                            conditionCurr = subject(sub).sess(sess).gem(gem).s(outcomeState).black_indices;
                        end
                        
                        for betas = 1:length(conditionCurr)

                            uniqueNameSpec = [uniqueName,num2str(conditionCurr(betas)),'*bf'];
                            index =[];
                            index = find(contains(SPM.xX.name,uniqueNameSpec)); %find indices of betas corresponding to contrast of interest
                            if isempty(index)
                            else
                            currentBetas = spm_read_vols(spm_vol(sprintf([currentModelFolder,...
                                fs, subPrefix,num2str(subjectString(sub)),fs,fileBase,sprintf('%03d',index),fileEnding]))); %read in beta file
                             
                            %put the relevant betas into a vector
                            currentBetas = currentBetas(:);
                            currentBetas = currentBetas(relevantVoxels.ROI); %select only those voxels that are part of your ROI
                            currentBetas = currentBetas(~isnan(currentBetas)); %in case there are edge-NaNs, exclude
                            currentVoxelsLocal(betas,:) = currentBetas'; %transpose (optional) and save in format currentVoxels(run,gem,betaValues) -> betas correspond to run
                            end 
                        end
                         currentVoxelsGlobal = [currentVoxelsGlobal;currentVoxelsLocal];
                    end
                    
                    meanCurrentVoxels = mean(currentVoxelsGlobal,1);
                    currentVoxelsReshaped((gem-1)*4+(outcomeState-1)*2 + door,:) = meanCurrentVoxels; %save as 16xvoxels matrix to be able to use pdist to obtain a 16x16 correlation matrix
                
                end
            end
            
        end
        
   
    % going to use the correlation distance (1-correlation coefficient) to
    % compute dissimilarity
    %for this, use reshaped currentVoxels matrix (16xvoxels matrix) -> pdist gives us
    %half of the 16x16 matrix of the pairwise dissimilarities (correlation distance) and squareform completes
    %the other triangle
    
    dissimLargeDep = squareform(pdist(currentVoxelsReshaped(1:8,:),'Correlation'));
    dissimLargeDepSub(sub,:,:) = dissimLargeDep(1:4,5:8); 
    simLargeDepSub(sub,:,:) = 1-dissimLargeDep(1:4,5:8); 
    
    dissimLargeIndep = squareform(pdist(currentVoxelsReshaped(9:16,:),'Correlation'));
    dissimLargeIndepSub(sub,:,:) = dissimLargeIndep(1:4,5:8);
    simLargeIndepSub(sub,:,:) = 1-dissimLargeIndep(1:4,5:8); 

end

%% now we will transform our data for plotting

for sub = 1:29
   %Calculate the mean zscored similarity on the diagonal (same choice &
   %outcome, different contexts) and off-diagonal (different choice & outcome,
   %different contexts); dep = dependent; indep = independent

    meanDiagDepZScored(sub) = squeeze(nanmean([atanh(1-dissimLargeDepSub(sub,1,1)),atanh(1-dissimLargeDepSub(sub,2,2)), atanh(1-dissimLargeDepSub(sub,3,3)),atanh(1-dissimLargeDepSub(sub,4,4))]));
    meanOffDiagDepZScored(sub) = squeeze(nanmean([atanh(1-dissimLargeDepSub(sub,1,2)),atanh(1-dissimLargeDepSub(sub,1,3)), atanh(1-dissimLargeDepSub(sub,1,4)),atanh(1-dissimLargeDepSub(sub,2,3)), ...
        atanh(1-dissimLargeDepSub(sub,2,4)),atanh(1-dissimLargeDepSub(sub,3,4)),atanh(1-dissimLargeDepSub(sub,2,1)),atanh(1-dissimLargeDepSub(sub,3,1)),atanh(1-dissimLargeDepSub(sub,3,2)),atanh(1-dissimLargeDepSub(sub,4,1)),atanh(1-dissimLargeDepSub(sub,4,2)),atanh(1-dissimLargeDepSub(4,3))]));

    meanDiagIndepZScored(sub) = squeeze(nanmean([atanh(1-dissimLargeIndepSub(sub,1,1)),atanh(1-dissimLargeIndepSub(sub,2,2)), atanh(1-dissimLargeIndepSub(sub,3,3)),atanh(1-dissimLargeIndepSub(sub,4,4))]));
    meanOffDiagIndepZScored(sub) = squeeze(nanmean([atanh(1-dissimLargeIndepSub(sub,1,2)),atanh(1-dissimLargeIndepSub(sub,1,3)), atanh(1-dissimLargeIndepSub(sub,1,4)),atanh(1-dissimLargeIndepSub(sub,2,3)), ...
        atanh(1-dissimLargeIndepSub(sub,2,4)),atanh(1-dissimLargeIndepSub(sub,3,4)),atanh(1-dissimLargeIndepSub(sub,2,1)),atanh(1-dissimLargeIndepSub(sub,3,1)),atanh(1-dissimLargeIndepSub(sub,3,2)),atanh(1-dissimLargeIndepSub(sub,4,1)),atanh(1-dissimLargeIndepSub(sub,4,2)),atanh(1-dissimLargeIndepSub(4,3))]));

end

%% boxplots with scatter (Figure 5)

dataset1 = [meanDiagDepZScored-meanOffDiagDepZScored];
dataset2 = [meanDiagIndepZScored-meanOffDiagIndepZScored];

colour = [.5,.5,.5];
colour2 = [225,225,225]/255;
xx = {dataset1; dataset2};
hold all;
boxplot([dataset1';dataset2'],[repmat(1,[29,1]);repmat(2,[29,1])],'Labels',{'Dependent', 'Independent'})
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colour,'FaceAlpha',.5);
end
xticks([1,2])
xticklabels({'Dependent','Independent'})
xtickangle(20)
scatterColours = [[24,116,205]/250,[165,42,42]/250];
for subs = 1:length(dataset1)
    plot(datasets,[dataset1(subs),dataset2(subs)], 'Color', colour2);
end
for k1 = 1:numel(datasets)
    scatter(ones(1,numel(xx{k1}))*datasets(k1), xx{k1}, 'filled','MarkerFaceColor',scatterColours((k1-1)*3+1:k1*3))
end
ylabel({'Difference in similarity between';'diagonal and off-diagonal'}, 'fontweight', 'bold');
set(gca,'FontSize',30)
set(gca,'XColor','black','YColor','black')
box 'off'

hold off

%% t-tests

%dependent condition
[h,pDep,ci,statsDep]= ttest(meanDiagDepZScored,meanOffDiagDepZScored);

%independent condition
[h,pIndep,ci,statsIndep]= ttest(meanDiagIndepZScored,meanOffDiagIndepZScored);

%interaction
[h,pInteraction,ci,statsInteraction] = ttest(meanDiagDepZScored-meanOffDiagDepZScored,meanDiagIndepZScored-meanOffDiagIndepZScored);

save(fileName);

end
