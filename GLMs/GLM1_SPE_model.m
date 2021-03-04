%% This script runs GLM1 - a GLM with two onset events (door and feedback)
%                          and four parametric modulators at the time of
%                          feedback: state prediction error, outcome
%                          valence, their interaction, and whether the
%                          current trial is a free or forced choice trial
%


% Loads into spm 1st level design matrix and estimates job
%
% cc: Steve Fleming 2008

% edited by Neil Garrett (NG) 08/19
% edited by Leonie Glitz (LG) 08/19
% edited by LG & NG 05/2020

clear all, clc;
cwd = pwd;

%load in onset file 
RLValues = load('/Volumes/Samsung_T5/gems/Outcome_Redo_May2020/RW_flex_onset.mat');


% Options for processing
compute_conditions = 1;
construct_design = 1;
estimate = 1;

spmdir = '/Users/leonieglitz/Documents/MATLAB/spm12/';

%% Change according to your directory strucutre and scan parameters
current_model_folder = 'Outcome_Redo_May2020/singleOutcomeOnset/';
dir_base    = '/Volumes/Samsung_T5/gems/'; %this is related to the above - work out how to balance paths between the two
fmri_folder = 'data/fMRI/processed/'; %this is where the processed fMRI files are 
dir_epi     = '/functional'; %this is where the epis are within each subject folder
sess_prfx   = '/session'; %this is the prefix for/name of the session folders in the dir_epi

num_subj       = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24,25,26,27,29,30,31];

fs = filesep;
blockNo = 4;
nslices = 32;
TR = 2.0000;
hpcutoff = 128; % hp filter

for n = 1:length(num_subj)
  
     
    subRL = RLValues.onset.sub(num_subj(n)); 
    
    % Go to subject directory    
     model_folder = [dir_base current_model_folder fs 'sub' num2str(num_subj(n))];
    
    % if the intended stats directory does not exits, make it and go into it
    if exist(model_folder) ~= 7
        mkdir(model_folder);
        cd(model_folder);
    end
    
    outputDir = model_folder;
    
    %=========================================================================
    %% Get onsets in scans from behavioural datafiles and define images
    %======================================================================
    
    for k = 1:blockNo;
        if compute_conditions
            disp(['Computing event onsets for subject ' num2str(num_subj(n)) ', session ' num2str(k)]);
            
            %subselect current block
            blockRL = subRL.sess(k);
            
            %calculate interaction of SPE and outcome
            interactionTerm =[]; %clear from previous block
            interactionTerm = blockRL.parametric.outcome .* blockRL.RLmod.SPE_abs;
            
            
            % Define behavioural data path
            risky_outcome =[];
            safe_outcome =[];
            names = [];
            onsets = [];
            durations = [];
            pmod = [];
            
            % three conditions, first control and the last two of interest
            names{1} = 'door_onset';
            onsets{1} = blockRL.categorical.door;
            durations{1} = 0;
            
           %feedback onsets when participants either gained or lost money
            names{2} = 'feedback';
            onsets{2} = blockRL.categorical.outcome;
            durations{2} =0;
            pmod(2).name{1} = 'outcome';
            pmod(2).param{1} = blockRL.parametric.outcome;
            pmod(2).poly{1} =1;
            pmod(2).name{2} = 'SPE';
            pmod(2).param{2} = blockRL.RLmod.SPE_abs;
            pmod(2).poly{2} =1;
            pmod(2).name{3} = 'interaction';
            pmod(2).param{3} = interactionTerm;
            pmod(2).poly{3} =1;
            pmod(2).name{4} = 'forcedTrial';
            pmod(2).param{4} = blockRL.parametric.force;
            pmod(2).poly{4} =1;
            
            % Use this to switch off orthogonalisation of pmods (only relevant if you have more than one per condition):
            orth{1} = 0;
            orth{2} =0 ;
            
            cd(outputDir);
            conditionFile = sprintf('conditions%d.mat', k);
            
            save(conditionFile, 'onsets', 'names', 'durations', 'pmod', 'orth');
            
        end
        
        %==========================================================================
        %% Construct design matrix
        %==========================================================================
        
        % Load files we have just created
        epiDir = [dir_base fmri_folder 'sub' num2str(num_subj(n)) dir_epi sess_prfx num2str(k)];

        mCorrFile = spm_select('List', epiDir, '^rp.*\.txt$');
        conditionFile = sprintf('conditions%d.mat',k);
        
        % Assign .mat file with onsets/names/pmods in to path
        conditionPath = [conditionFile];
        multiregPath = [epiDir fs mCorrFile];
        
        % select scans and concatenate
        f   = spm_select('List', epiDir, '^swraf.*\.img$');     % Select smoothed normalised images
        files  = cellstr([repmat([epiDir fs],size(f,1),1) f]);
        % clear temporary variables for next run
        jobs{1}.stats{1}.fmri_spec.sess(k).scans = files;
        f = []; files = [];
        
        jobs{1}.stats{1}.fmri_spec.sess(k).multi = {conditionPath};
        jobs{1}.stats{1}.fmri_spec.sess(k).multi_reg = {multiregPath};
        % high pass filter
        jobs{1}.stats{1}.fmri_spec.sess(k).hpf = hpcutoff;
        jobs{1}.stats{1}.fmri_spec.sess(k).cond = struct([]);
        jobs{1}.stats{1}.fmri_spec.sess(k).regress = struct([]);
    end
    
    %==========================================================================
    %======================================================================
    jobs{1}.stats{1}.fmri_spec.dir = {outputDir};
    % timing variables
    jobs{1}.stats{1}.fmri_spec.timing.units     = 'secs';
    jobs{1}.stats{1}.fmri_spec.timing.RT        = TR;
    jobs{1}.stats{1}.fmri_spec.timing.fmri_t    = nslices;
    jobs{1}.stats{1}.fmri_spec.timing.fmri_t0   = nslices/2;
    
    % basis functions
    jobs{1}.stats{1}.fmri_spec.bases.hrf.derivs         = [0 0]; %%added the derivatives of the HRF
    % model interactions (Volterra) OPTIONS: 1|2 = order of convolution
    jobs{1}.stats{1}.fmri_spec.volt                     = 1;
    % global normalisation
    jobs{1}.stats{1}.fmri_spec.global                   = 'None';
    % explicit masking
    jobs{1}.stats{1}.fmri_spec.mask                     = {[spmdir fs 'tpm/mask_ICV.nii']};
    % serial correlations
    jobs{1}.stats{1}.fmri_spec.cvi                      = 'AR(1)';
    % no factorial design
    jobs{1}.stats{1}.fmri_spec.fact = struct('name', {}, 'levels', {});
    
    
    %==========================================================================
    %% run model specification
    %==========================================================================
    if construct_design
%         cd(outputDir);
        % save and run job
        save specify.mat jobs
        disp(['RUNNING model specification for subject ' num2str(num_subj(n))]);
        spm_jobman('run','specify.mat');
        clear jobs
    end
    
    % Ensure implicit masking is switched off
    load SPM
    SPM.xM.TH = repmat(-Inf, length(SPM.xM.TH), 1);
    SPM.xM.I = 0;
    save SPM SPM
    
    %% Estimate
    % setup job structure for model estimation and estimate
    % ---------------------------------------------------------------------
    if estimate
        jobs{1}.stats{1}.fmri_est.spmmat = {[outputDir fs 'SPM.mat']};
        save estimate.mat jobs
        disp(['RUNNING model estimation for subject ' num2str(num_subj(n))])
        spm_jobman('run','estimate.mat');
        clear jobs
    end
end
cd(cwd);
