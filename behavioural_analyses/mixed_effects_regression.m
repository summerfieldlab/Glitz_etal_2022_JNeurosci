%% This script will run a 5-back lagged regression analysis on the behavioural data.
%More specifically,it will use information from the previous five trials coded in terms of
%whether participants should choose the black door on the current trial
%given previous trials to predict choice behaviour. If participants had
%understood the block structure, we would expect information from the other
%Context to be more predictive of choice behaviour in dependent than in
%independent blocks. 


%in a second step, it will then average the information content from the
%same and the other context across the previous five trials and take the difference between the two 
%to assess whether the average information value significantly interacts with block type in
%predicting choice behaviour on choice trials

clc;
clear all;
close all;

data = readtable('gem_dat.csv');

%remove exclusions
data = data(data.participantID~=21,:);
data = data(data.participantID~=28,:);

%preallocate empty arrays
nBack_participant = [];
nBackSameContext_participant = [];
nBackOtherContext_participant = [];
nBack =[];
nBackSameContext =[];
nBackOtherContext =[];

%loop over all participants
for participant = 1:max(data.participantID)
    
    %pull data for that participant 
    if participant == 21 | participant == 28 %exclude those excluded from fMRI analysis
    else 
        %clear data from previous participant
        nBack_participant = [];
        nBackSameContext_participant = [];
        nBackOtherContext_participant = [];
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
                    
                    %sort into same and other context as presented on current trial
                    if currentData.gem_presented(trial) == currentData.gem_presented(trial-trialsBack)
                        nBackSameContext_block(trial,trialsBack) = nBack_block(trial,trialsBack);
                        nBackOtherContext_block(trial,trialsBack) = 0; %other Context not seen
                    
                    else %not the same Context trialsBack trials before as on current trial
                        nBackSameContext_block(trial,trialsBack) = 0; %same context not seen
                        nBackOtherContext_block(trial,trialsBack) = nBack_block(trial,trialsBack);
                    end 
                end   
            end
            
            %concatenate blocks together for current participant
            nBack_participant = [nBack_participant;nBack_block];
            nBackSameContext_participant = [nBackSameContext_participant;nBackSameContext_block];
            nBackOtherContext_participant = [nBackOtherContext_participant;nBackOtherContext_block];
            
        end 
        
        %concatenate data together across participants
        nBack = [nBack; nBack_participant];
        nBackSameContext = [nBackSameContext; nBackSameContext_participant];
        nBackOtherContext = [nBackOtherContext; nBackOtherContext_participant];
        
    end
end

%% does information from the same and other context on the previous five trials predict choices (on choice trials)?

%convert pick_black to 1s and -1s for regression
chooseBlack = data.pick_black;
% chooseBlack(chooseBlack == 0) = -1; 
participantID = data.participantID(data.participantID~=21 & data.participantID~=28);
blockType = data.blockType(data.participantID~=21 & data.participantID~=28);
blockType(blockType == 2) = -1; %independent blocks as -1 rather than 2
forcedTrial = data.forcedTrial;

%make the predictors into proper variables
oneBackSameContext = nBackSameContext(:,1);
twoBackSameContext = nBackSameContext(:,2);
threeBackSameContext = nBackSameContext(:,3);
fourBackSameContext = nBackSameContext(:,4);
fiveBackSameContext = nBackSameContext(:,5);

oneBackOtherContext = nBackOtherContext(:,1);
twoBackOtherContext = nBackOtherContext(:,2);
threeBackOtherContext = nBackOtherContext(:,3);
fourBackOtherContext = nBackOtherContext(:,4);
fiveBackOtherContext = nBackOtherContext(:,5);

oneBack = nBack(:,1);
twoBack = nBack(:,2);
threeBack = nBack(:,3);
fourBack = nBack(:,4);
fiveBack = nBack(:,5);


%save relevant values in a table so that fitglm can work with them
GLMTable = table(participantID,chooseBlack,oneBack, twoBack,threeBack,fourBack,...
    fiveBack,oneBackSameContext,twoBackSameContext,threeBackSameContext,fourBackSameContext,...
    fiveBackSameContext, oneBackOtherContext,twoBackOtherContext,threeBackOtherContext,fourBackOtherContext,...
    fiveBackOtherContext,blockType, forcedTrial);

%% 
%run multiple logistic regression analysis for all block types
allBlocks = fitglme(GLMTable(GLMTable.forcedTrial==0,:), ...
    'chooseBlack ~ oneBack + twoBack + threeBack + fourBack + fiveBack + (1 + oneBack + twoBack + threeBack + fourBack+ fiveBack | participantID)',...
    'Distribution','binomial');

%split data into dependent and independent blocks

%run multiple logistic regression analysis for dependent blocks
dependentBlocks = fitglme(GLMTable(GLMTable.forcedTrial==0 & GLMTable.blockType == 1,:), ...
    'chooseBlack ~ oneBackSameContext + twoBackSameContext + threeBackSameContext + fourBackSameContext + fiveBackSameContext + oneBackOtherContext + twoBackOtherContext + threeBackOtherContext + fourBackOtherContext + fiveBackOtherContext + (1 +oneBackSameContext + twoBackSameContext + threeBackSameContext + fourBackSameContext + fiveBackSameContext + oneBackOtherContext + twoBackOtherContext + threeBackOtherContext + fourBackOtherContext + fiveBackOtherContext | participantID)',...
    'Distribution','binomial');

%run multiple logistic regression analysis for independent blocks
independentBlocks = fitglme(GLMTable(GLMTable.forcedTrial==0 & GLMTable.blockType == -1,:), ...
    'chooseBlack ~ oneBackSameContext + twoBackSameContext + threeBackSameContext + fourBackSameContext + fiveBackSameContext + oneBackOtherContext + twoBackOtherContext + threeBackOtherContext + fourBackOtherContext + fiveBackOtherContext + (1 +oneBackSameContext + twoBackSameContext + threeBackSameContext + fourBackSameContext + fiveBackSameContext + oneBackOtherContext + twoBackOtherContext + threeBackOtherContext + fourBackOtherContext + fiveBackOtherContext | participantID)',...
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

