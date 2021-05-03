function [pDep, statsDep,ciDep, pIndep, statsIndep,ciIndep, pInteraction, statsInteraction, ciInteraction] = RSA_outcome(mask, fileName)
%% This function is going to run an RSA at the time of outcome in the ROI 
%  that is specified. It will split trials up by the door chosen (light/dark) 
%  and whether participants reached the gain/loss or the safe state. Then 
%  we will compare the representational similarity of identical and 
%  non-identical action-outcome pairings across the two contexts in 
%  each of our conditions (dependent and independent). 

%  cc: Leonie Glitz, University of Oxford, 2020

%  (a similar script for a predecessor of this analysis was written by Neil
%  Garrett too)

% TAKES INPUTS: mask (i.e. mask name) and fileName -> save all data
% under what name 

%% Step 1: load the masked voxel data

load(['/Users/leonieglitz/Desktop/Garrett&Glitz_etal_2021/outcome_fMRI_data/BOLD_outcome_',mask,'.mat']);
%this loads a struct
%BOLD(sub).sess(sess).gem(context).state(outcomeState).white_indices or .black_indices
%where data is a trial x voxel matrix with the outputs of GLM2a for the
%mask of interest, sorted by subject, session, context (==gem) and
%probability range (subjective probabilities of presented door transitioning
%to riksy state, divided into quartiles)

%% Step 2: compute the representational (dis)similarity between contexts and door chosen - outcome state combinations

%loop through subjects
for sub = 1:29
    currentVoxelsReshaped = []; 
    currentVoxelsGlobal =[];
        
    %loop through conditions/gems and find corresponding beta.nii files
        for context = 1:4 %conditions
            for outcomeState = 1:2
                for door = 1:2
                    
                    currentVoxelsGlobal =[]; 
                    for sess = 1:4
                        
                        BOLD_data = [];                       
                        if door == 1
                            BOLD_data = BOLD(sub).sess(sess).gem(context).state(outcomeState).light_indices;
                        else
                            BOLD_data = BOLD(sub).sess(sess).gem(context).state(outcomeState).dark_indices;
                        end
                        
                        currentVoxelsGlobal = [currentVoxelsGlobal;BOLD_data];
                    end
                    
                    meanCurrentVoxels = mean(currentVoxelsGlobal,1);
                    currentVoxelsReshaped((context-1)*4+(outcomeState-1)*2 + door,:) = meanCurrentVoxels; %save as 16xvoxels matrix to be able to use pdist to obtain a 16x16 correlation matrix
                
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
datasets = [1,2];

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
[h,pDep,ciDep,statsDep]= ttest(meanDiagDepZScored,meanOffDiagDepZScored);

%independent condition
[h,pIndep,ciIndep,statsIndep]= ttest(meanDiagIndepZScored,meanOffDiagIndepZScored);

%interaction
[h,pInteraction,ciInteraction,statsInteraction] = ttest(meanDiagDepZScored-meanOffDiagDepZScored,meanDiagIndepZScored-meanOffDiagIndepZScored);

save(fileName);

end
