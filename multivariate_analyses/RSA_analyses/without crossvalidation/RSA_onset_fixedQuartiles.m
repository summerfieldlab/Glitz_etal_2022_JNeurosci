%% This function is going to run an RSA at the time of door onset in the ROI
%  that is specified. It will split trials up into four probability bins by their
%  probability of the presented door leading to the gain/loss state. Then
%  we will compare the representational similarity of identical and
%  non-identical transition probability bins across the two contexts in
%  each of our conditions (dependent and independent).

%  cc: Leonie Glitz, University of Oxford, 2020

%please note: this script reads in the original data at the moment, so
%cannot be executed using the data provided online (we cannot make this
%data available due to ethics restrictions)

%% Step 1: calculate beta indices of conditions of interest for each sub!

load('/Users/leonieglitz/Downloads/RW_flex_onset.mat');
behDat = readtable('/Volumes/Samsung_T5/gems/gem_dat.csv');

behDatNoNaN = behDat(~isnan(behDat.chooseLeft),:);

includedSubs = [1:20 22:27 29:31];


frequencyCount(1:length(includedSubs),1:16,1) = 0;


for sub = 1:length(includedSubs)
    currSubject = onset.sub(includedSubs(sub));
    
    for sess = 1:4
        
        currSession = currSubject.sess(sess);
        
        for gemNo = 1:4
            
            %make below into structs so that diff trial nums do not matter
            
            gem_dat_Data = find(currSession.RLmod.prob_rl_doorpresented_combined <= .25 & currSession.parametric.gem_n==gemNo);
            frequencyCount(sub,(4*(gemNo-1)+1)) = frequencyCount(sub,(4*(gemNo-1)+1)) + length(gem_dat_Data);
            subject(sub).sess(sess).gem(gemNo).prob(1).data = gem_dat_Data;
            
            gem_dat_Data = [];
            
            gem_dat_Data = find(currSession.RLmod.prob_rl_doorpresented_combined > .25 & currSession.RLmod.prob_rl_doorpresented_combined <= .5 & currSession.parametric.gem_n==gemNo);
            frequencyCount(sub,(4*(gemNo-1)+2)) = frequencyCount(sub,(4*(gemNo-1)+2)) + length(gem_dat_Data);
            subject(sub).sess(sess).gem(gemNo).prob(2).data = gem_dat_Data;
            
            gem_dat_Data = [];
            
            gem_dat_Data = find(currSession.RLmod.prob_rl_doorpresented_combined > .5 & currSession.RLmod.prob_rl_doorpresented_combined <= .75 & currSession.parametric.gem_n==gemNo);
            frequencyCount(sub,(4*(gemNo-1)+3)) = frequencyCount(sub,(4*(gemNo-1)+3)) + length(gem_dat_Data);
            subject(sub).sess(sess).gem(gemNo).prob(3).data = gem_dat_Data;
            
            gem_dat_Data = [];
            
            gem_dat_Data = find(currSession.RLmod.prob_rl_doorpresented_combined > .75 & currSession.RLmod.prob_rl_doorpresented_combined <= 1 & currSession.parametric.gem_n==gemNo);
            frequencyCount(sub,(4*(gemNo-1)+4)) = frequencyCount(sub,(4*(gemNo-1)+4)) + length(gem_dat_Data);
            subject(sub).sess(sess).gem(gemNo).prob(4).data = gem_dat_Data;
            
            gem_dat_Data = [];
        end
        
        
    end
    
    
    
end


%% Step 2: read in the betas corresponding to the different conditions and perform RSA

fs = filesep; %so script can be used on mac or windows

%setup
currentModelFolder = '/Volumes/Samsung_T5/gems/singleTrialModelNeil/Onset'; %where your models are generally

subPrefix = 'sub'; %note: to go through all subjects, you need to take the spmT_000X.nii file from the individual subject folders
% as there are several files with similar names, it is important to cd into
% it rather than the 2nd level results

numRuns = 4;

subjectString = [1:20 22:27 29:31]; %this is important if you have exclusions

fileBase = 'beta_0'; %base for reading in the beta images
fileEnding = '.nii'; %ending
subPrefix = 'sub';

currentVoxelsReshaped = [];

relevantVoxels = struct(); %ensures that you don't get a struct error later
relevantVoxels.mask = spm_read_vols(spm_vol('SPE_mask_neg.nii')); %read in mask for ROI
relevantVoxels.ROI = find(relevantVoxels.mask); %find where ROI is 1 in the mask - these are the indices of the betas from the  ROI we want

%loop through subjects
for sub = 1:length(subjectString)
    currentVoxelsReshaped =[];
    currentVoxelsGlobal =[];
    %load the SPM.mat file
    load([currentModelFolder,...
        fs,subPrefix,num2str(subjectString(sub)),fs,'SPM.mat']); %load in the current subject's SPM file
    
    %loop through gems and find corresponding beta.nii files
    for gem = 1:4 %here: blocks
        for probRange = 1:4
            currentBetas = [];
            currentVoxels = [];
            currentVoxelsLocal =[];
            currentVoxelsGlobal =[];
            for sess = 1:4
                uniqueName = ['Sn(',num2str(sess),') door_onset_'];
                currentVoxelsLocal = [];
                conditionCurr = subject(sub).sess(sess).gem(gem).prob(probRange).data; %here sub is used for indexing as data encoded in struct above indexed by sub
                
                for betas = 1:length(conditionCurr)
                    
                    uniqueNameSpec = [uniqueName,num2str(conditionCurr(betas)),'*bf'];
                    index =[];
                    index = find(contains(SPM.xX.name,uniqueNameSpec)); %find index of beta corresponding to contrast of interest
                    if isempty(index)
                    else
                        currentBetas = spm_read_vols(spm_vol(sprintf([currentModelFolder,...
                            fs, subPrefix,num2str(subjectString(sub)),fs,fileBase,sprintf('%03d',index),fileEnding]))); %read in beta file
                        
                        %put the relevant betas into a vector
                        currentBetas = currentBetas(:);
                        currentBetas = currentBetas(relevantVoxels.ROI); %select only those voxels that are part of your ROI
                        currentBetas = currentBetas(~isnan(currentBetas)); %in case there are edge-NaNs, exclude
                        currentVoxelsLocal(betas,:) = currentBetas'; %transpose (optional) and save in format currentVoxels(trial,betaValues)
                    end
                end
                currentVoxelsGlobal = [currentVoxelsGlobal;currentVoxelsLocal]; %concatenate across sessions/blocks
            end
            
            meanCurrentVoxels = nanmean(currentVoxelsGlobal,1); %mean across trials
            
            if isempty(meanCurrentVoxels)
                'empty mean condition' %if this happens, think about inserting NaN values 
                currentVoxelsReshaped((gem-1)*4+probRange,:) = repmat(NaN,[1,69]);
            else
                currentVoxelsReshaped((gem-1)*4+probRange,:) = meanCurrentVoxels; %save as 16xvoxels matrix to be able to use pdist to obtain a 16x16 correlation matrix
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
    simLargeDepSub(sub,:,:) = 1- dissimLargeDep(1:4,5:8);
    
    dissimLargeIndep = squareform(pdist(currentVoxelsReshaped(9:16,:),'Correlation'));
    dissimLargeIndepSub(sub,:,:) = dissimLargeIndep(1:4,5:8);
    simLargeIndepSub(sub,:,:) = 1-dissimLargeIndep(1:4,5:8);
end

%% now we will transform our data for plotting

for sub = 1:29
    
    %Calculate the mean zscored similarity on the diagonal (same prob.
    %bins, different contexts) and off-diagonal (different probability
    %bins, different contexts); dep = dependent; indep = independent
    
    meanDiagDepZScored(sub) = squeeze(nanmean([atanh(1-dissimLargeDepSub(sub,1,1)),atanh(1-dissimLargeDepSub(sub,2,2)), atanh(1-dissimLargeDepSub(sub,3,3)),atanh(1-dissimLargeDepSub(sub,4,4))]));
    meanOffDiagDepZScored(sub) = squeeze(nanmean([atanh(1-dissimLargeDepSub(sub,1,2)),atanh(1-dissimLargeDepSub(sub,1,3)), atanh(1-dissimLargeDepSub(sub,1,4)),atanh(1-dissimLargeDepSub(sub,2,3)), ...
        atanh(1-dissimLargeDepSub(sub,2,4)),atanh(1-dissimLargeDepSub(sub,3,4)),atanh(1-dissimLargeDepSub(sub,2,1)),atanh(1-dissimLargeDepSub(sub,3,1)),atanh(1-dissimLargeDepSub(sub,3,2)),atanh(1-dissimLargeDepSub(sub,4,1)),atanh(1-dissimLargeDepSub(sub,4,2)),atanh(1-dissimLargeDepSub(4,3))]));
    
    meanDiagIndepZScored(sub) = squeeze(nanmean([atanh(1-dissimLargeIndepSub(sub,1,1)),atanh(1-dissimLargeIndepSub(sub,2,2)), atanh(1-dissimLargeIndepSub(sub,3,3)),atanh(1-dissimLargeIndepSub(sub,4,4))]));
    meanOffDiagIndepZScored(sub) = squeeze(nanmean([atanh(1-dissimLargeIndepSub(sub,1,2)),atanh(1-dissimLargeIndepSub(sub,1,3)), atanh(1-dissimLargeIndepSub(sub,1,4)),atanh(1-dissimLargeIndepSub(sub,2,3)), ...
        atanh(1-dissimLargeIndepSub(sub,2,4)),atanh(1-dissimLargeIndepSub(sub,3,4)),atanh(1-dissimLargeIndepSub(sub,2,1)),atanh(1-dissimLargeIndepSub(sub,3,1)),atanh(1-dissimLargeIndepSub(sub,3,2)),atanh(1-dissimLargeIndepSub(sub,4,1)),atanh(1-dissimLargeIndepSub(sub,4,2)),atanh(1-dissimLargeIndepSub(4,3))]));
    
    
end
%% boxplots with scatter (Figure 3 & 4)

dataset1 = [meanDiagDepZScored-meanOffDiagDepZScored];
dataset2 = [meanDiagIndepZScored-meanOffDiagIndepZScored];


colour = [.5,.5,.5];
colour2 = [225,225,225]/255;

hold all;
boxplot([dataset1';dataset2'],[repmat(1,[29,1]);repmat(2,[29,1])],'Labels',{'Dependent', 'Independent'})
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colour,'FaceAlpha',.5);
end

datasets = [1,2];
xx = {dataset1; dataset2};
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


%% run stats on the data

%difference between diagonal and off diagonal similarity in the dependent
%condition
[h,pDep,ciDep,statsDep]= ttest(meanDiagDepZScored,meanOffDiagDepZScored)

%difference between diagonal and off-diagonal similarity in the independent
%condition
[h,pIndep,ciIndep,statsIndep]= ttest(meanDiagIndepZScored,meanOffDiagIndepZScored);

%interaction
[h,pInteraction,ciInteraction,statsInteraction] = ttest(meanDiagDepZScored-meanOffDiagDepZScored,meanDiagIndepZScored-meanOffDiagIndepZScored)
