%% This script is going to run a whole-brain searchlight RSA at the time of outcome presentation
%  that is specified. It will split trials up by the door chosen (light/dark) 
%  and whether participants reached the gain/loss or the safe state. Then 
%  we will compare the representational similarity of identical and 
%  non-identical action-outcome pairings across the two contexts in 
%  each of our conditions (dependent and independent). 

%  cc: Leonie Glitz, University of Oxford, 2020

%  NOTE: this script uses functions from and is based on functions in the RSA toolbox,
%        which can be found here: http://www.mrc-cbu.cam.ac.uk/methods-and-resources/toolboxes/license/
%        Nili et al. (2014). A toolbox for representational similarity analysis. 
%        PLoS computational biology, 10(4), e1003553.

%  for this script to run, the RSA toolbox folder needs to be added to the
%  working directory
%% Step 1: calculate beta indices of conditions of interest for each sub!

load('/Users/leonieglitz/Downloads/RW_flex_onset.mat');
behDat = readtable('/Volumes/Samsung_T5/gems/gem_dat.csv');

behDatNoNaN = behDat(~isnan(behDat.chooseLeft),:);

includedSubs = [1:20 22:27 29:31];

%we want to run separate RSAs for dependent and independent contexts with
%themselves; within that we want to split up by door colour (here denoted
%black and white) and outcome state (s)

%find indices of relevant trials
for sub = 1:29
    currSubject = behDatNoNaN(behDatNoNaN.participantID == includedSubs(sub),:);
    for sess = 1:4
        currSession = currSubject(currSubject.block_n == sess+1,:);
        
        %indices for context 1
        context_dat_Data = find(currSession.pick_black==0 & currSession.gem_presented==1 & currSession.outcomeState==1);
        subject(sub).sess(sess).context(1).s(1).white_indices = context_dat_Data;
        
        context_dat_Data = find(currSession.pick_black==1 & currSession.gem_presented==1 & currSession.outcome~=0);
        subject(sub).sess(sess).context(1).s(1).black_indices = context_dat_Data;
        
        context_dat_Data = find(currSession.pick_black~=1 & currSession.gem_presented==1 & currSession.outcome==0);
        subject(sub).sess(sess).context(1).s(2).white_indices = context_dat_Data;
        
        context_dat_Data = find(currSession.pick_black==1 & currSession.gem_presented==1 & currSession.outcome==0);
        subject(sub).sess(sess).context(1).s(2).black_indices = context_dat_Data;
        
        
        %indices for context 2
        context_dat_Data = find(currSession.pick_black~=1 & currSession.gem_presented==2 & currSession.outcome~=0);
        subject(sub).sess(sess).context(2).s(1).white_indices = context_dat_Data;

        context_dat_Data = find(currSession.pick_black==1 & currSession.gem_presented==2 & currSession.outcome~=0);
        subject(sub).sess(sess).context(2).s(1).black_indices = context_dat_Data;
        
        context_dat_Data = find(currSession.pick_black~=1 & currSession.gem_presented==2 & currSession.outcome==0);
        subject(sub).sess(sess).context(2).s(2).white_indices = context_dat_Data;
        
        context_dat_Data = find(currSession.pick_black==1 & currSession.gem_presented==2 & currSession.outcome==0);
        subject(sub).sess(sess).context(2).s(2).black_indices = context_dat_Data;
        
        
        %indices for context 3
        context_dat_Data = find(currSession.pick_black~=1 & currSession.gem_presented==3 & currSession.outcome~=0);
        subject(sub).sess(sess).context(3).s(1).white_indices = context_dat_Data;
        
        context_dat_Data = find(currSession.pick_black==1  & currSession.gem_presented==3 & currSession.outcome~=0);
        subject(sub).sess(sess).context(3).s(1).black_indices = context_dat_Data;
        
        context_dat_Data = find(currSession.pick_black~=1  & currSession.gem_presented==3 & currSession.outcome==0);
        subject(sub).sess(sess).context(3).s(2).white_indices = context_dat_Data;
        
        context_dat_Data = find(currSession.pick_black==1 & currSession.gem_presented==3 & currSession.outcome==0);
        subject(sub).sess(sess).context(3).s(2).black_indices = context_dat_Data;
        
        
        %indices for context 4
        context_dat_Data = find(currSession.pick_black~=1 & currSession.gem_presented==4 & currSession.outcome~=0)
        subject(sub).sess(sess).context(4).s(1).white_indices = context_dat_Data;
        
        context_dat_Data = find(currSession.pick_black==1 & currSession.gem_presented==4 & currSession.outcome~=0);
        subject(sub).sess(sess).context(4).s(1).black_indices = context_dat_Data;
        
        context_dat_Data = find(currSession.pick_black~=1 & currSession.gem_presented==4 & currSession.outcome==0);
        subject(sub).sess(sess).context(4).s(2).white_indices = context_dat_Data;
        
        context_dat_Data = find(currSession.pick_black==1 & currSession.gem_presented==4 & currSession.outcome==0);
        subject(sub).sess(sess).context(4).s(2).black_indices = context_dat_Data;
    end
    
 
    
end

%% Step 2: format data so that borrowed toolbox scripts can read in the data (for the one modified RSA toolbox script, see toolbox functions subfolder)

%SETUP FOR TOOLBOX FUNCTION

relevantVoxels = struct(); %ensures that you don't get a struct error later
relevantVoxels.mask = spm_read_vols(spm_vol('/Volumes/Samsung_T5/gems/RSA/gems_EC/wholeBrainAttempt.img')); %read in the mask for each sub
relevantVoxels.ROI = find(relevantVoxels.mask); %find where ROI is 1 in the mask - these are the betas from the  ROI we want

mask = spm_read_vols(spm_vol('/Volumes/Samsung_T5/gems/RSA/gems_EC/wholeBrainAttempt.img')); %read in the mask again for toolbox function later
maskName = 'wholeBrain';

%model RDMs 
load('/Volumes/Samsung_T5/gems/singleTrialModel/outcomeOnlyModel.mat');
outcomeOnlyDep = outcomeOnlyModel; 
models(1).RDM = outcomeOnlyDep;
models(1).name = 'outcomeDependent';

load('/Volumes/Samsung_T5/gems/singleTrialModel/outcomeOnlyModelIndep.mat');
indepOnly = outcomeOnlyModel;
models(2).RDM  = indepOnly;
models(2).name = 'outcomeIndependent'; 

%control model RDMs looking for greater similarity for the same action only
load('/Volumes/Samsung_T5/gems/singleTrialModelNeil/Outcome/actionOnlyDep.mat')
models(3).RDM = actionOnlyDep;
models(3).name = 'sameActionDependent';

load('/Volumes/Samsung_T5/gems/singleTrialModelNeil/Outcome/actionOnlyIndep.mat')
models(4).RDM = actionOnlyIndep;
models(4).name = 'sameActionIndependent';

%control model RDMs looking for greater similarity for the same outcome
%state only
load('/Volumes/Samsung_T5/gems/singleTrialModelNeil/Outcome/stateOnlyDep.mat')
models(5).RDM = stateOnlyDep;
models(5).name = 'sameStateDependent';

load('/Volumes/Samsung_T5/gems/singleTrialModelNeil/Outcome/stateOnlyIndep.mat')
models(6).RDM = stateOnlyIndep;
models(6).name = 'sameStateIndependent';


%import the RSA toolbox functions we need
import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*


userOptions.analysisName = 'outcomeSearchlight';
userOptions.rootPath = '/Users/leonieglitz/Desktop/Garrett&Glitz_etal_2021/multivariate_analyses/RSA_analyses/RSA_searchlights/outcomeSearchlight/';
userOptions.maskNames = {'whole_brain'};
userOptions.voxelSize = [3.5 3.5 3.5];
userOptions.searchlightRadius = 10.5;
userOptions.averageSessions = 0; 
localOptions.monitor = 0; 
localOptions.averageSessions = 0; 
localOptions.fisher = 0; 


%SETUP FOR READING IN BETAS

fs = filesep; %so script can be used on mac or windows

%setup
currentModelFolder = '/Volumes/Samsung_T5/gems/singleTrialModelNeil/Outcome'; %where your models are generally

subPrefix = 'sub'; %note: to go through all subjects, you need to take the spmT_000X.nii file from the individual subject folders
% as there are several files with similar names, it is important to cd into
% it rather than the 2nd level results

numRuns = 4;

subjectString = [1:20 22:27 29:31]; %this is important if you have exclusions

fileBase = 'beta_0'; %base for reading in the beta images
fileEnding = '.nii'; %ending
subPrefix = 'sub';


currentVoxelsReshaped = [];

%read in the betas for each subject
for sub = 1:length(subjectString)
    
    currentVoxelsReshaped =[];
    currentVoxelsGlobal =[];
    meanCurrentVoxels =[];
    
    %load the SPM.mat file
    load([currentModelFolder,...
        fs,subPrefix,num2str(subjectString(sub)),fs,'SPM.mat']); %load in the current subject's SPM file
    
    subjectMetadataStruct = spm_vol([currentModelFolder,...
        fs,subPrefix,num2str(subjectString(sub)),fs,'beta_0001.nii']);
    
    %loop through contexts and find corresponding beta.nii files
        for context = 1:4 
            for outcomeState = 1:2
                for door = 1:2

                    currentBetas = [];
                   
                    currentVoxelsLocal =[];
                    currentVoxelsGlobal =[]; 
                    
                    for sess = 1:4
                        currentVoxelsLocal = [];
                        uniqueName = ['Sn(',num2str(sess),') outcome_onset_'];
                        
                        if door == 1
                            conditionCurr = subject(sub).sess(sess).context(context).s(outcomeState).white_indices;
                        else
                            conditionCurr = subject(sub).sess(sess).context(context).s(outcomeState).black_indices;
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
                            %currentBetas = currentBetas(~isnan(currentBetas)); %in case there are edge-NaNs, exclude
                            currentVoxelsLocal(betas,:) = currentBetas'; %transpose (optional) and save in format currentVoxels(run,context,betaValues) -> betas correspond to run
                            end 
                        end
                         currentVoxelsGlobal = [currentVoxelsGlobal;currentVoxelsLocal];
                    end
                    
                    meanCurrentVoxels = mean(currentVoxelsGlobal,1);
                    currentVoxelsReshaped(:,(context-1)*4+(outcomeState-1)*2 + door,1) = meanCurrentVoxels; %save as voxels x 16 matrix to be able to feed into the searchlightMapping_fMRI function

                     
                end
            end
            
        end
        
        %feed the voxels into the searchlight function (see toolbox
        %functions subfolder)
        subject4Searchlight = ['sub',num2str(subjectString(sub))]; 
   
        [rs, ps, ns, searchlightRDMs.(subject4Searchlight)]  = searchlightMapping_fMRI(currentVoxelsReshaped, models, mask, userOptions, localOptions);
        
        nMaps_nS = struct();
        rMaps_nS = struct(); 
        nMaps_nS.(subject4Searchlight).(maskName) = ns(:,:,:); % How many voxels contributed to the searchlight centred at each point. (Those with n==1 are excluded because the results aren't multivariate.)

    for modelNumber = 1:numel(models)

        modelName = (models(modelNumber).name);
        
        % Store results in indexed volumes
        rMaps_nS.(modelName).(subject4Searchlight).(maskName) = rs(:,:,:,modelNumber); % r-values for correlation with each model RDM

        %% Save native space version

        % Write the native-space r-map to a file
        rMapMetadataStruct_nS = subjectMetadataStruct;
        rMapMetadataStruct_nS.fname = fullfile(userOptions.rootPath, '/','Maps', '/',strcat(userOptions.analysisName, '_rMap_', maskName, '_', modelName, '_', subject4Searchlight, '.img'));
        rMapMetadataStruct_nS.descrip =  'R-map';
        rMapMetadataStruct_nS.dim = size(rMaps_nS.(modelName).(subject4Searchlight).(maskName));

        cd([userOptions.rootPath, 'Maps']);

        rsa.spm.spm_write_vol(rMapMetadataStruct_nS, rMaps_nS.(modelName).(subject4Searchlight).(maskName));
    end

    clear fullBrainVolumes rs ps ns;

    %% Save relevant info
    % The analysisName will be used to label the files which are eventually saved.
    mapsFilename = [userOptions.analysisName, '_fMRISearchlight_Maps.mat'];
    RDMsFilename = [userOptions.analysisName, '_fMRISearchlight_RDMs.mat'];
    DetailsFilename = [userOptions.analysisName, '_fMRISearchlight_Details.mat'];
    timeStamp = datestr(now);

    fprintf(['Saving searchlight maps']);
    save([userOptions.rootPath,'Maps',mapsFilename], 'rMaps_nS', 'nMaps_nS', '-v7.3');

    fprintf(['Saving RDMs to ' fullfile(userOptions.rootPath, 'RDMs', RDMsFilename) '\n']);
    save([userOptions.rootPath,'RDMs',RDMsFilename], 'searchlightRDMs', '-v7.3');

    fprintf(['Saving Details to ' fullfile(userOptions.rootPath, 'Details', DetailsFilename) '\n']);
    save([userOptions.rootPath,'Details',DetailsFilename], 'timeStamp', 'userOptions', '-v7.3');

    printProgress = ['finished subject ', num2str(includedSubs(sub))]
end
