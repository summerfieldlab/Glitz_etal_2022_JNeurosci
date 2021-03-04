%% This script runs GLM2a - it creates single-trial onset timepoint events 
%                           to be used for the subsequent RSA analysis 

%  template cc: Steve Fleming 2008 (metalab)
%  edited by NG / LG 08/19
%  edited by Neil Garrett, 06/2020
%  edited by Leonie Glitz, 06/2020 

%script to get design matrix onsets
%Loads into spm 1st level design matrix and estimates job

%clear everything
clear all; close all; clc;

%store present working dir
cwd = pwd;

%load onset file
load('RW_flex_onset.mat');

%options for processing
compute_conditions = 1;
construct_design = 1;
estimate = 1;

%ensure spm directory is in the path
spmdir = '/Users/leonieglitz/Documents/MATLAB/spm12';
addpath(spmdir);

%% Change according to your directory strucutre and scan parameters
dir_base    = '/Volumes/Samsung_T5/gems';
curr_model_path = '/singleTrialModelNeil/Onset';
scans_folder = 'data/fMRI/processed'; 
dir_epi     = 'functional';
sess_prfx   = 'session';

sub_ns       = [1:20 22:27 29:31];

fs = filesep;

%number of functional runs
blockNo = 4;

nslices = 32;
TR = 2;
hpcutoff = 128; %hp filter

for n = 1:length(sub_ns)
 
     %go to subject directory
    model_folder = [dir_base fs curr_model_path fs 'sub' num2str(sub_ns(n))];
    
    %if the intended stats directory does not exits, make it and cd to it
    if (exist(model_folder) ~= 7)
        mkdir(model_folder);
        cd(model_folder);
    else
    %if it does exist, go into it and delete existing contrast files
    %as we are re-estimating
        cd(model_folder);
        delete('SPM.mat', '*.img', '*.hdr');
    end
    
    outputDir = model_folder;

    for k = 1:blockNo;
        
        if compute_conditions

        %=========================================================================
        %% Get onsets in scans from behavioural datafiles and define images
        %=========================================================================
        
            disp(['Computing event onsets for subject ' num2str(sub_ns(n)) ', session ' num2str(k)]);
            names       = []; %important to delete everything!!!
            onsets      = [];
            durations   = [];
            pmod        = [];
            orth = [];
            
            %door onset screen - loop over
            for dOns = 1:length(onset.sub(sub_ns(n)).sess(k).categorical.door)
                
                names{dOns} = ['door_onset_',num2str(dOns)];
                onsets{dOns} = onset.sub(sub_ns(n)).sess(k).categorical.door(dOns);
                durations{dOns} = 0; %stick function
                orth{dOns} = 0; %off
                
            end 
            

                names{length(onset.sub(sub_ns(n)).sess(k).categorical.door)+1} = 'outcome_onset';
                onsets{length(onset.sub(sub_ns(n)).sess(k).categorical.door)+1} = onset.sub(sub_ns(n)).sess(k).categorical.outcome;
                durations{length(onset.sub(sub_ns(n)).sess(k).categorical.door)+1} = 0; %delta (stick) function
                orth{length(onset.sub(sub_ns(n)).sess(k).categorical.door)+1} = 0; %off



            cd(outputDir);
            conditionFile = sprintf('conditions%d.mat', k);
            
            save(conditionFile, 'onsets', 'names', 'durations', 'pmod', 'orth');
            
        end
        
        %==========================================================================
        %% Construct design matrix
        %==========================================================================
        
        %path to where functional scans are...
        epiDir = [dir_base fs scans_folder fs 'sub' num2str(sub_ns(n)) fs dir_epi fs sess_prfx num2str(k)];
                
        %pull out movement parameters (to add to design mtrix as a nuisence)
        mCorrFile = spm_select('List', epiDir, '^rp.*\.txt$');
        conditionFile = sprintf('conditions%d.mat', k);
        
        %assign .mat file with onsets/names/pmods in to path
        conditionPath = [outputDir fs conditionFile];
        multiregPath = [epiDir fs mCorrFile];
        
        %select scans and concatenate
        f = spm_select('List', epiDir, '^swraf.*\.img$');     % Select smoothed normalised images
        files  = cellstr([repmat([epiDir fs], size(f, 1), 1) f]);
        
        %clear temporary variables for next run
        jobs{1}.stats{1}.fmri_spec.sess(k).scans = files;
        f = []; files = [];
        
        jobs{1}.stats{1}.fmri_spec.sess(k).multi = {conditionPath};
        jobs{1}.stats{1}.fmri_spec.sess(k).multi_reg = {multiregPath};
        
        %high pass filter
        jobs{1}.stats{1}.fmri_spec.sess(k).hpf = hpcutoff;
        jobs{1}.stats{1}.fmri_spec.sess(k).cond = struct([]);
        jobs{1}.stats{1}.fmri_spec.sess(k).regress = struct([]);
        
    end
    
    jobs{1}.stats{1}.fmri_spec.dir = {outputDir};
    
    %timing variables
    jobs{1}.stats{1}.fmri_spec.timing.units     = 'secs';
    jobs{1}.stats{1}.fmri_spec.timing.RT        = TR;
    jobs{1}.stats{1}.fmri_spec.timing.fmri_t    = nslices;
    jobs{1}.stats{1}.fmri_spec.timing.fmri_t0   = nslices/2;
    
    %basis functions
    jobs{1}.stats{1}.fmri_spec.bases.hrf.derivs         = [0 0];
    
    %model interactions (Volterra) OPTIONS: 1|2 = order of convolution
    jobs{1}.stats{1}.fmri_spec.volt                     = 1;
    
    %global normalisation
    jobs{1}.stats{1}.fmri_spec.global = 'None';
    
    %explicit masking
    jobs{1}.stats{1}.fmri_spec.mask = {[spmdir fs 'tpm/mask_ICV.nii']};
    
    %serial correlations
    jobs{1}.stats{1}.fmri_spec.cvi = 'AR(1)';
    
    %no factorial design
    jobs{1}.stats{1}.fmri_spec.fact = struct('name', {}, 'levels', {});
    
    %==========================================================================
    %% run model specification
    %==========================================================================
    if construct_design
        
        cd(outputDir);
        
        %save and run job
        save specify.mat jobs
        disp(['RUNNING model specification for subject ' num2str(sub_ns(n))]);
        spm_jobman('run','specify.mat');
        clear jobs
        
    end
    
    %ensure implicit masking is switched off
    load SPM
    SPM.xM.TH = repmat(-Inf, length(SPM.xM.TH), 1);
    SPM.xM.I = 0;
    save SPM SPM
    
    %% Estimate
    
    %setup job structure for model estimation and estimate
    if estimate
        jobs{1}.stats{1}.fmri_est.spmmat = {[outputDir fs 'SPM.mat']};
        save estimate.mat jobs
        disp(['RUNNING model estimation for subject ' num2str(sub_ns(n))]);
        spm_jobman('run','estimate.mat');
        clear jobs
    end
    
end

cd(cwd);
