%% This script is going to run a whole-brain searchlight RSA at the time of door onset 
%  that is specified. It will split trials up into quartiles by their 
%  probability of the presented door leading to the gain/loss state. Then 
%  we will compare the representational similarity of identical and 
%  non-identical transition probability bins across the two contexts in 
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
%themselves; within that we want to split up by model-based prob of
%gain/loss of door presented


frequencyCount(1:length(includedSubs),1:16,1) = 0; 


for sub = 1:length(includedSubs)
    currSubject = onset.sub(includedSubs(sub));
    
    currProbTmp = [currSubject.sess(1).RLmod.prob_rl_doorpresented_combined;currSubject.sess(2).RLmod.prob_rl_doorpresented_combined;currSubject.sess(3).RLmod.prob_rl_doorpresented_combined;currSubject.sess(4).RLmod.prob_rl_doorpresented_combined];
    currProbEsts = sort(currProbTmp);
    medianVal = median(currProbEsts);
    quartilesCurr = [currProbEsts(round(length(currProbEsts)/4)),medianVal,currProbEsts(round(3* length(currProbEsts)/4)) , max(currProbEsts)]; 
    
    for sess = 1:4
        
        currSession = currSubject.sess(sess);
                
        for contextNo = 1:4
        
        %make below into structs so that diff trial nums do not matter
        
        indicesTmp = find(currSession.RLmod.prob_rl_doorpresented_combined <= quartilesCurr(1) & currSession.parametric.gem_n==contextNo);
        frequencyCount(sub,(4*(contextNo-1)+1)) = frequencyCount(sub,(4*(contextNo-1)+1)) + length(indicesTmp); 
        subject(sub).sess(sess).context(contextNo).prob(1).data = indicesTmp;
        
        indicesTmp = [];

        indicesTmp = find(currSession.RLmod.prob_rl_doorpresented_combined > quartilesCurr(1) & currSession.RLmod.prob_rl_doorpresented_combined <= quartilesCurr(2) & currSession.parametric.gem_n==contextNo);
        frequencyCount(sub,(4*(contextNo-1)+2)) = frequencyCount(sub,(4*(contextNo-1)+2)) + length(indicesTmp); 
        subject(sub).sess(sess).context(contextNo).prob(2).data = indicesTmp;
        
        indicesTmp = [];
        
        indicesTmp = find(currSession.RLmod.prob_rl_doorpresented_combined > quartilesCurr(2) & currSession.RLmod.prob_rl_doorpresented_combined <= quartilesCurr(3) & currSession.parametric.gem_n==contextNo);
        frequencyCount(sub,(4*(contextNo-1)+3)) = frequencyCount(sub,(4*(contextNo-1)+3)) + length(indicesTmp); 
        subject(sub).sess(sess).context(contextNo).prob(3).data = indicesTmp;
        
        indicesTmp = [];
        
        indicesTmp = find(currSession.RLmod.prob_rl_doorpresented_combined > quartilesCurr(3) & currSession.RLmod.prob_rl_doorpresented_combined <= 1 & currSession.parametric.gem_n==contextNo);
        frequencyCount(sub,(4*(contextNo-1)+4)) = frequencyCount(sub,(4*(contextNo-1)+4)) + length(indicesTmp); 
        subject(sub).sess(sess).context(contextNo).prob(4).data = indicesTmp;
        
        indicesTmp = [];
        end 

       
    end
    
    
end


%% Step 2: format data so that borrowed toolbox scripts can read in the data (for the one modified RSA toolbox script, see toolbox functions subfolder)

%SETUP FOR TOOLBOX FUNCTION

relevantVoxels = struct(); %ensures that you don't get a struct error later
relevantVoxels.mask = spm_read_vols(spm_vol('/Volumes/Samsung_T5/gems/RSA/gems_EC/wholeBrainAttempt.img')); %read in the mask for each sub
relevantVoxels.ROI = find(relevantVoxels.mask); %find where ROI is 1 in the mask - these are the betas from the  ROI we want

%read in the mask used - here: whole-brain
mask = spm_read_vols(spm_vol('/Volumes/Samsung_T5/gems/RSA/gems_EC/wholeBrainAttempt.img')); %read in the mask for each sub
maskName = 'wholeBrain';

%read in the model RDMs
load('/Users/leonieglitz/Desktop/Garrett&Glitz_etal_2021/multivariate_analyses/RSA_analyses/RSA_searchlights/model RDMs/outcomeOnlyModel.mat');
outcomeOnlyDep = outcomeOnlyModel; 
models(1).RDM = outcomeOnlyDep;
models(1).name = 'outcomeDependent';

load('/Users/leonieglitz/Desktop/Garrett&Glitz_etal_2021/multivariate_analyses/RSA_analyses/RSA_searchlights/model RDMs/outcomeOnlyModelIndep.mat');
indepOnly = outcomeOnlyModel;
models(2).RDM  = indepOnly;
models(2).name = 'outcomeIndependent'; 


%import the toolbox functions we need
import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

%set toolbox settings
userOptions.analysisName = 'onsetSearchlight';
userOptions.rootPath = '/Users/leonieglitz/Desktop/Garrett&Glitz_etal_2021/multivariate_analyses/RSA_analyses/RSA_searchlights/onsetSearchlight/';
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
currentModelFolder = '/Volumes/Samsung_T5/gems/singleTrialModelNeil/Onset'; %where your models are generally

numRuns = 4;

subjectString = [1:20 22:27 29:31]; %subjects to include

fileBase = 'beta_0'; %base for reading in the beta images
fileEnding = '.nii'; %ending
subPrefix = 'sub';


currentVoxelsReshaped = [];

%loop through subjects
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
            for probRange = 1:4
                    
                    
                    currentBetas = [];
                    currentVoxels = [];
                    currentVoxelsLocal =[];
                    currentVoxelsGlobal =[]; 
                    for sess = 1:4
                        uniqueName = ['Sn(',num2str(sess),') door_onset_'];
                        currentVoxelsLocal = []; 
                        conditionCurr = subject(sub).sess(sess).context(context).prob(probRange).data;
                        
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
                            currentVoxelsLocal(betas,:) = currentBetas'; %transpose (optional) and save in format currentVoxels(run,context,betaValues) -> betas correspond to run
                            end 
                        end
                         currentVoxelsGlobal = [currentVoxelsGlobal;currentVoxelsLocal];
                    end
                    
                    NaNArray(1:43955,1) = NaN;
            
                    if isempty(currentVoxelsGlobal)
                        'currentVoxelsGlobal is empty'
                        currentVoxelsReshaped(:,(context-1)*4+probRange,1) = NaNArray; 
                    else
                        meanCurrentVoxels = mean(currentVoxelsGlobal,1);
                        currentVoxelsReshaped(:,(context-1)*4+probRange,1) = meanCurrentVoxels; %save as voxels x 16 matrix to be able to feed into the searchlightMapping_fMRI function
                    end 
                     
             end
            
        end
        
        %feed the voxel values into the searchlight function
        subject4Searchlight = ['sub',num2str(subjectString(sub))]; 
   
        [rs, ps, ns, searchlightRDMs.(subject4Searchlight)]  = searchlightMapping_fMRI(currentVoxelsReshaped, models, mask, userOptions, localOptions);
        
        nMaps_nS = struct();
        rMaps_nS = struct(); 
        nMaps_nS.(subject4Searchlight).(maskName) = ns(:,:,:); % How many voxels contributed to the searchlight centred at each point. (Those with n==1 are excluded because the results aren't multivariate.)

    for modelNumber = 1:numel(models)

        modelName = (models(modelNumber).name);
        
        % Store results in indexed volumes
        rMaps_nS.(modelName).(subject4Searchlight).(maskName) = rs(:,:,:,modelNumber); % spearman r-values for correlation with each model RDM

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
