% -------------------------------------------
%    Author: Serafeim Loukas, July 13, 2021
%    For PPI-CAPs project with Joana (& Lorena)
%
%    Updated & modified by: Serafeim Loukas on July 20, 2021
% -------------------------------------------
%
% NOTE: If paths are not correctly defined in SPM.mat use "spm_change_paths" function to change them.
%
% Expected output 
%  ------------------------------------------------------------------------
% 21-Jul-2021 14:03:50 - Running job #2
% ------------------------------------------------------------------------
% 21-Jul-2021 14:03:50 - Running 'Volume of Interest'
%    VOI saved as /Users/loukas/Desktop/PPI_CAPs/subjects/PTMT1_49RiGo_Rechercheneonat_3151/results_mask_new/Voice_fMRI_2019_RPARAM24/VOI_Auditory_1.mat
% 21-Jul-2021 14:03:50 - Done    'Volume of Interest'
% 21-Jul-2021 14:03:51 - Done

clc;
clear all;

%% General paths setup
% set main path (universal path)
% -------------------------------
ppicapsDir = '/Users/loukas/Desktop/PPI_CAPs/';
%ppicapsDir = '/home/loukas/Project_PPI_CAPs/NeoPPI-CAPs/';

% where to save the CAPs
save_caps = [ppicapsDir 'results/'];
if not(isfolder(save_caps)); mkdir(save_caps);end

% Go to subject folder
addpath(genpath('/home/loukas/spm12/'));
addpath(genpath([ppicapsDir, 'scripts/']));
mainPath = [ppicapsDir, 'subjects/'];
cd(mainPath)

S = dir('*T1*');
dirFlags = [S.isdir];
subFolders = S(dirFlags);

%% Create cell with the subject folder names
% -------------------------------
n_subjects= size(S,1);
subject_name= {};
for k = 1 : n_subjects
    subject_name{k} = subFolders(k).name;
end

%% VOI creation
 for k = 1 : n_subjects
     % path to subjects
     path_to_subj = [mainPath subject_name{k} filesep];
     path_to_spm = [path_to_subj 'results_mask_new' filesep 'Voice_fMRI_2019_RPARAM24' filesep];
     dir_spm = strvcat(path_to_spm);
     path_to_motion = [path_to_subj 'functional'];
     
     % get and load the SPM.mat file
     [spm_file] = spm_select('FPList',deblank(path_to_spm),['SPM.mat']);
     
     % select correct file and avoid weird ._ files on server
     if size(spm_file,1) > 1
         [d,b,c]    = fileparts(spm_file(1,:));
         [d2,b2,c2] = fileparts(spm_file(2,:));
         if length(b)==3
             spm_file = deblank(spm_file(1,:));
         elseif length(b2)==3
             spm_file = deblank(spm_file(2,:));      
         end
     end
     
     spm_f=load(spm_file);
     nSessions = length(spm_f.SPM.nscan);
     
     for thisSession = 1:nSessions         
         matlabbatch{1}.spm.util.voi.spmmat = cellstr(fullfile(deblank(spm_file)));
         matlabbatch{1}.spm.util.voi.adjust = 1;
         matlabbatch{1}.spm.util.voi.session = thisSession;
         matlabbatch{1}.spm.util.voi.name = 'Auditory';
         % Auditory3 is in T2 space but will automatically resliced to functional space by SPM
         matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {[ppicapsDir 'subjects/Auditory3.nii,1']}; 
         % Auditory3 is in the template T2 neonatal space. SPM will convert that and it 
         % will create the final "VOI_Auditory_1.mat" with 499 voxels in the functional space.
         % see VOI_Auditory_mask.nii (that is in functional space with 499 nnz voxels)
         matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
         matlabbatch{1}.spm.util.voi.expression = 'i1';
         spm_jobman('run',matlabbatch);
     end
    
 end
 % NOTE: If paths are not correctly defined in SPM.mat use "spm_change_paths" function to change them.
 % Better alternative is to create the VOI_Auditory_1.mat once and use its coordinates for all subjects.
 % Data are registered to the same space so the 2 ways are equivalent. 
    
