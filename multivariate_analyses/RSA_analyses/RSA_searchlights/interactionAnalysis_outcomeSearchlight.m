%% This script is going to read in the neural RDMs for the dependent and independent condition at time of outcome
%  and calculate an interaction - i.e. where in the brain the difference
%  between on- and off-diagonal similarity is greater in the dependent than
%  the independent condition. It will read in the RDMs saved when running
%  outcomeSearchlight.m

%  c: Leonie Glitz, University of Oxford, 2020


analysisName = 'outcomeSearchlight';



%don't make into triangular matrix for Spearmans correlation because the
%neural RDM is not symmetric around the diagonal!

load(['/Volumes/Samsung_T5/gems/singleTrialModelNeil/Outcome/RDMs',analysisName,'_fMRISearchlight_RDMs.mat']); 
includedSubs = [1:20 22:27 29:31];

for sub = 1:numel(fieldnames(searchlightRDMs))
    currentSub = searchlightRDMs.(['sub',num2str(includedSubs(sub))]);
    for voxelX = 1:45
        for voxelY = 1:55
            for voxelZ = 1:45
                
                %get neural RDM for current voxel
                currentNeuralRDM = squeeze(currentSub(:,:,voxelX,voxelY,voxelZ));
                
                %do all computations on the Fisher's z transformed
                %similarity scores (sim == 1-dissim)
                
                %calculate difference diag-offDiag for dependent
                %between-gems -> NOTE: include both sides of the off-diagonal as these are NOT symmetric!
                depDiagMinOffDiag(sub,voxelX,voxelY,voxelZ) = nanmean([atanh(1-currentNeuralRDM(5,1)),atanh(1-currentNeuralRDM(6,2)),atanh(1-currentNeuralRDM(7,3)),atanh(1-currentNeuralRDM(8,4))])... %diagonal mean
                    - nanmean([atanh(1-currentNeuralRDM(5,2)),atanh(1-currentNeuralRDM(5,3)), atanh(1-currentNeuralRDM(5,4)),atanh(1-currentNeuralRDM(6,1)),atanh(1-currentNeuralRDM(6,3)),...
                    atanh(1-currentNeuralRDM(6,4)),atanh(1-currentNeuralRDM(7,1)),atanh(1-currentNeuralRDM(7,2)), atanh(1-currentNeuralRDM(7,4)),atanh(1-currentNeuralRDM(8,1)),atanh(1-currentNeuralRDM(8,2)),atanh(1-currentNeuralRDM(8,3))]);
         

                %calculate difference diag-offDiag for independent between-gems -> NOTE: include both sides of the off-diagonal as these are NOT symmetric!       
                indepDiagMinOffDiag(sub,voxelX,voxelY,voxelZ) = nanmean([atanh(1-currentNeuralRDM(9,13)),atanh(1-currentNeuralRDM(10,14)),atanh(1-currentNeuralRDM(11,15)),atanh(1-currentNeuralRDM(12,16))])... %diagonal mean
                    - nanmean([atanh(1-currentNeuralRDM(9,14)),atanh(1-currentNeuralRDM(9,15)), atanh(1-currentNeuralRDM(9,16)),atanh(1-currentNeuralRDM(10,13)),atanh(1-currentNeuralRDM(10,15)), ...
                    atanh(1-currentNeuralRDM(10,16)),atanh(1-currentNeuralRDM(11,13)),atanh(1-currentNeuralRDM(11,14)),atanh(1-currentNeuralRDM(11,16)),atanh(1-currentNeuralRDM(12,13)),atanh(1-currentNeuralRDM(12,14)),atanh(1-currentNeuralRDM(12,15))]);

                %save difference maps
                differenceMap(sub,voxelX,voxelY,voxelZ) = depDiagMinOffDiag(sub,voxelX,voxelY,voxelZ) - indepDiagMinOffDiag(sub,voxelX,voxelY,voxelZ);                
            
            end
        end
    end 
    
    currentDiffMeans = squeeze(differenceMap(sub,:,:,:)); 
    
    info = spm_vol('/Volumes/Samsung_T5/gems/RSA/visual_disc_gainLoss_doorColour/sub24/beta_0001.nii');
    info = rmfield(info,'descrip');
    info = rmfield(info,'private');
    info = rmfield(info,'fname');
    
    newInfo.dt = info.dt;
    newInfo.mat = info.mat;
    newInfo.dim = size(currentDiffMeans);
    newInfo.pinfo = info.pinfo;
    newInfo.fname = fullfile(['/Volumes/Samsung_T5/gems/singleTrialModelNeil/Outcome/Interaction/differenceMap_',analysisName,'_sub',num2str(sub),'.nii']);
    newInfo.descrip = 'tmap';
    

    V = spm_write_vol(newInfo,currentDiffMeans);

    
    ['Finished writing difference maps for subject ',num2str(sub),'/',num2str(length(includedSubs))]
    
end

