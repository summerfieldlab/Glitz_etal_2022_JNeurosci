%% This script will run a 5-back lagged regression analysis on the behavioural data of the fMRI version of the heist task.
%More specifically,it will use information from the previous five trials coded in terms of
%whether participants should choose the black door on the current trial
%given previous trials to predict choice behaviour. If participants had
%understood the block structure, we would expect information from the other
%gem to be more predictive of choice behaviour in dependent than in
%independent blocks. 


%in a second step, it will then average the information content from the
%same and the other gem across the previous five trials and take the difference between the two 
%to assess whether the average information value significantly interacts with block type in
%predicting choice behaviour on choice trials

clc;
clear all;
close all;

data = readtable('data/gem_dat_fMRI.csv');

%remove exclusions
% data = data(~isnan(data.chooseLeft),:);
data = data(data.participantID~=21,:);
data = data(data.participantID~=28,:);

%preallocate empty arrays
nBack_participant = [];
nBackSameGem_participant = [];
nBackOtherGem_participant = [];
nBack =[];
nBackSameGem =[];
nBackOtherGem =[];

%loop over all participants
for participant = 1:max(data.participantID)
    
    %pull data for that participant 
    if participant == 21 | participant == 28 %exclude those excluded from fMRI analysis
    else 
        %clear data from previous participant
        nBack_participant = [];
        nBackSameGem_participant = [];
        nBackOtherGem_participant = [];
        partData = data(data.participantID == participant,:);
        
        %pull data for current block
        for blocks = 1:max(data.block_n)
            currentData = partData(partData.block_n == blocks,:);
            %% main loop where regressors are created
            for trial = 6:height(currentData)
                %loop over the previous five trials- assume you want to get
                %to gain/loss state for now and then reverse later if loss trial
                for trialsBack = 1:5
                    if currentData.pick_black(trial-trialsBack) == 1 %chose dark door
                        %outcome of choosing black
                        if currentData.outcome(trial-trialsBack) == 0 %safe state
                            nBack_block(trial,trialsBack) = -1; %the trialsBackth trial back does not suggest black leads to gainloss
                        else %risky state
                            nBack_block(trial,trialsBack) = 1; %trialsBackth trial does suggest black leads to gain/loss
                        end 
                    elseif isnan(currentData.pick_black(trial-trialsBack))
                        nBack_block(trial,trialsBack) = NaN;
                    else %chose white door
                        %outcome of choosing light door
                        if currentData.outcome(trial-trialsBack) == 0 %safe state
                            nBack_block(trial,trialsBack) = 1; %the trialsBackth trial suggests black leads to gainloss (as white does not)
                        else %risky state
                            nBack_block(trial,trialsBack) = -1; %trialsBackth trial does not suggest black leads to gain/loss (because white does)
                        end 
                    end
                    
                    if currentData.rew_loss(trial) == 1 %gain trial
                        %keep sign - no change
                    else %loss trial
                        %reverse sign - do not want to reach gain/loss state
                        nBack_block(trial,trialsBack) = nBack_block(trial,trialsBack)*(-1);
                    end
                    
                    %sort into same and other gem as presented on current trial
                    if currentData.gem_presented(trial) == currentData.gem_presented(trial-trialsBack)
                        nBackSameGem_block(trial,trialsBack) = nBack_block(trial,trialsBack);
                        nBackOtherGem_block(trial,trialsBack) = 0; %other gem not seen
                    else %not the same gem trialsBack trials before as on current trial
                        nBackSameGem_block(trial,trialsBack) = 0; %same gem not seen
                        nBackOtherGem_block(trial,trialsBack) = nBack_block(trial,trialsBack);
                    end 
                end   
            end
            
            %concatenate blocks together for current participant
            nBack_participant = [nBack_participant;nBack_block];
            nBackSameGem_participant = [nBackSameGem_participant;nBackSameGem_block];
            nBackOtherGem_participant = [nBackOtherGem_participant;nBackOtherGem_block];
            
        end 
        
        %concatenate data together across participants
        nBack = [nBack; nBack_participant];
        nBackSameGem = [nBackSameGem; nBackSameGem_participant];
        nBackOtherGem = [nBackOtherGem; nBackOtherGem_participant];
        
    end
end

%% does information from the same and other gem on the previous five trials predict choices (on choice trials)?

%convert pick_black to 1s and -1s for regression
chooseBlack = data.pick_black;
% chooseBlack(chooseBlack == 0) = -1; 
participantID = data.participantID(data.participantID~=21 & data.participantID~=28);
blockType = data.blockType(data.participantID~=21 & data.participantID~=28);
blockType(blockType == 2) = -1; %independent blocks as -1 rather than 2
forcedTrial = data.forcedTrial;

%make the predictors into proper variables
oneBackSameGem = nBackSameGem(:,1);
twoBackSameGem = nBackSameGem(:,2);
threeBackSameGem = nBackSameGem(:,3);
fourBackSameGem = nBackSameGem(:,4);
fiveBackSameGem = nBackSameGem(:,5);

oneBackOtherGem = nBackOtherGem(:,1);
twoBackOtherGem = nBackOtherGem(:,2);
threeBackOtherGem = nBackOtherGem(:,3);
fourBackOtherGem = nBackOtherGem(:,4);
fiveBackOtherGem = nBackOtherGem(:,5);

oneBack = nBack(:,1);
twoBack = nBack(:,2);
threeBack = nBack(:,3);
fourBack = nBack(:,4);
fiveBack = nBack(:,5);


%save relevant values in a table so that fitglm can work with them
GLMTable = table(participantID,chooseBlack,oneBack, twoBack,threeBack,fourBack,...
    fiveBack,oneBackSameGem,twoBackSameGem,threeBackSameGem,fourBackSameGem,...
    fiveBackSameGem, oneBackOtherGem,twoBackOtherGem,threeBackOtherGem,fourBackOtherGem,...
    fiveBackOtherGem,blockType, forcedTrial);

%% 
%run multiple logistic regression analysis for all block types
allBlocks = fitglme(GLMTable(GLMTable.forcedTrial==0,:), ...
    'chooseBlack ~ oneBack + twoBack + threeBack + fourBack + fiveBack + (1 + oneBack + twoBack + threeBack + fourBack+ fiveBack | participantID)',...
    'Distribution','binomial');

%split data into dependent and independent blocks

%run multiple logistic regression analysis for dependent blocks
dependentBlocks = fitglme(GLMTable(GLMTable.forcedTrial==0 & GLMTable.blockType == 1,:), ...
    'chooseBlack ~ oneBackSameGem + twoBackSameGem + threeBackSameGem + fourBackSameGem + fiveBackSameGem + oneBackOtherGem + twoBackOtherGem + threeBackOtherGem + fourBackOtherGem + fiveBackOtherGem + (1 +oneBackSameGem + twoBackSameGem + threeBackSameGem + fourBackSameGem + fiveBackSameGem + oneBackOtherGem + twoBackOtherGem + threeBackOtherGem + fourBackOtherGem + fiveBackOtherGem | participantID)',...
    'Distribution','binomial');

%run multiple logistic regression analysis for independent blocks
independentBlocks = fitglme(GLMTable(GLMTable.forcedTrial==0 & GLMTable.blockType == -1,:), ...
    'chooseBlack ~ oneBackSameGem + twoBackSameGem + threeBackSameGem + fourBackSameGem + fiveBackSameGem + oneBackOtherGem + twoBackOtherGem + threeBackOtherGem + fourBackOtherGem + fiveBackOtherGem + (1 +oneBackSameGem + twoBackSameGem + threeBackSameGem + fourBackSameGem + fiveBackSameGem + oneBackOtherGem + twoBackOtherGem + threeBackOtherGem + fourBackOtherGem + fiveBackOtherGem | participantID)',...
    'Distribution','binomial');

%% plot the results 

%create a plot for dependent blocks including errorbars (SEM)
figure;
scatter([1:5]+0.1,dependentBlocks.Coefficients(2:6,2),'r','filled');
hold on
b = errorbar([1:5]+0.1,dependentBlocks.Coefficients(2:6,2),dependentBlocks.Coefficients.SE(2:6),'LineStyle','none');
b.Color ='black';
hold on
h = errorbar([1:5]-0.1,dependentBlocks.Coefficients(7:11,2),dependentBlocks.Coefficients.SE(7:11),'LineStyle','none');
h.Color ='black';
hold on
scatter([1:5]-0.1,dependentBlocks.Coefficients(7:11,2),'b','filled');
xlim([0,6]);
xlabel('Number of trials back');
ylabel('Regression coefficient (beta)');
title('Influence of CC (red) and NCVC (blue) Evidence on Choices in Dependent Blocks');
ylim([-0.4,1.4]);

%create a plot for independent blocks including errorbars (SEM)
figure;
scatter([1:5]+0.1,independentBlocks.Coefficients(2:6,2),'r','filled');
hold on
b = errorbar([1:5]+0.1,independentBlocks.Coefficients(2:6,2),independentBlocks.Coefficients.SE(2:6),'LineStyle','none');
b.Color ='black';
hold on
h = errorbar([1:5]-0.1,independentBlocks.Coefficients(7:11,2),independentBlocks.Coefficients.SE(7:11),'LineStyle','none');
h.Color ='black';
hold on
scatter([1:5]-0.1,independentBlocks.Coefficients(7:11,2),'b','filled');
xlim([0,6]);
xlabel('Number of trials back');
ylabel('Regression coefficient (beta)');
title('Influence of CC (red) and NCVC (blue) Evidence on Choices in Independent Blocks');
ylim([-0.4,1.4]);

%% 
%in a second step, we then average the information content from the
%same and the other gem across the previous five trials and take the differ
%ence between the two values to assess whether the average information value
%significantly interacts with block type in predicting choice behaviour on choice trials

%calculate the averages of the same and other gem information and take the
%difference between them
for inclTrials = 1:height(GLMTable)
    averageSame = mean(nBackSameGem(inclTrials,:),'omitnan');
    averageOther = mean(nBackOtherGem(inclTrials,:),'omitnan');
    evidenceDifference(inclTrials,1) = averageSame - averageOther;
end

%put that difference into a table with the other relevant bits of information
differenceTable = table(evidenceDifference, blockType, chooseBlack, participantID,forcedTrial);

%run a regression on that 
differenceRegression = fitglme(differenceTable(differenceTable.forcedTrial==0,:),'chooseBlack ~ evidenceDifference*blockType + (1 + evidenceDifference * blockType |participantID)','Distribution','binomial');

% save 'dependentBlocks.mat'
% writematrix(dependentBlocks.Coefficients.Estimate(2:11)', 'stats_dependent_estimates_data.csv');
% writematrix(dependentBlocks.Coefficients.SE(2:11)', 'stats_dependent_se_data.csv');
% 
% save 'independentBlocks.mat'
% writematrix(independentBlocks.Coefficients.Estimate(2:11)', 'stats_independent_estimates_data.csv');
% writematrix(independentBlocks.Coefficients.SE(2:11)', 'stats_independent_se_data.csv');

% anova(differenceRegression,'dfmethod','satterthwaite')
% 
% anova(differenceRegression,'dfmethod','residual')
% Satterthwaite
% 
% writetable(GLMTable, 'mixed_model_r_dat_lagged.csv')
% 
% writetable(differenceTable, 'mixed_model_r_dat_averaged.csv')
