% -------------------------------------------------------
%
%    main - function to execute PPI-CAPs workflow
%
%    Created by:         Lorena Freitas 
%    Last checked:       28.09.2019
%
%   
% ------------------------------------------------------

function main

% ------------------------------------------------------
% Set up
% ------------------------------------------------------

% Distance used for kmeans 
dist        = 'cosine';

% Subjects in analysis
sujets      = {'Sujet01', 'Sujet10', 'Sujet16', 'Sujet17', 'Sujet19' ...
               'Sujet21', 'Sujet26', 'Sujet27', 'Sujet31', 'Sujet32', ...
               'Sujet09', 'Sujet14', 'Sujet20', 'Sujet23', 'Sujet25' };

ppicapsDir = '<PATH TO SAVE PPI-CAPs>';  
 
% ------------------------------------------------------
% Load data for all subjects
% ------------------------------------------------------

% Initialise variables
supraFrames=[];subjInfo=[];timeInfo=[];taskInfo=[]; seedSignLabels=[]; 

% Loop through subjects' data to make super subject. The following
% information must have been saved for each subject: the suprathreshold
% frames (selected based on seed's most (de)active moments after whole brain
% deconvolution), the indeces of each suprathreshold frame within the 
% original time course (timeInfo), and task labels for each suprathreshold 
% frame indicating which task was being performed at that point in time 
% (taskInfo); 
for i=1:length(sujets)
    
    % Load subject's variables, paths, etc
    thisSubject = initialize_vars(sujets{i}, 'control');
    
    % Select suprathresholdFrames
    [PPIframes, suprathresholdFramesIx, correspondingTask4Frames] = thresholdFrames(thisSubject, 'PCC', 60);
    load([thisSubject.dataDir thisSubject.curSubj '.mat']); % loads PPIframes, suprathresholdFramesIx and correspondingTask4Frames for each subject
    
    supraFrames      = [supraFrames, PPIframes']; % suprathreshold frames from each subject in the shape #voxels x #frames
    subjInfo         = [subjInfo, ones(1,size( PPIframes',2))*i]; % subject indeces corresponding to each supraframe
    timeInfo         = [timeInfo, suprathresholdFramesIx]; % each frame's corresponding index within a subject's timecourse
    taskInfo         = [taskInfo, correspondingTask4Frames]; % task each frame belongs to
    
    seedSignLabels = [seedSignLabels suprathresholdFramesSeedSign];
    
end

% Find voxels in and out of Grey Matter
[voxels, nonVoxels] = findGMvoxels;

% Find matrix rotations for saving nifti files
[voxel_size, voxel_shift] = prepareToSaveNifti;


% ------------------------------------------------------
% Cluster timepoints, using frames as feature vectors
% ------------------------------------------------------

% Choose number of clusters (PPI-CAPs)
k = 4;

% Run kmeans++ on 
[PPI_CAPs, PPI_CAPs_Ix]  = kmeansCustom(supraFrames(voxels,:), k, 100); 
   
for i = 1:k
        CAPix{i}              = find(PPI_CAPs == i); % #frames x 1; index of frames for this CAP
        CAPl{i}               = PPI_CAPs == i;
        scf{i}                = CAPl{i}.* PPI_CAPs_Ix; %signed cluster's frames
        cluster{i}            = supraFrames(:,CAPix{i}); % #voxels x #frames
        CAPmean               = (sum(supraFrames(:,scf{i}==1),2) - sum(supraFrames(:,scf{i}==2),2))/sum(CAPl{i} );
        CAPmean(nonVoxels,:)  = 0;
        CAPmean{i}            = CAPmean;
        
        polarityLabels= scf{i};
        polarityLabels = polarityLabels(CAPix{i});
        taskLabels = taskInfo(CAPix{i});  
        polarityLabels{i} = polarityLabels;
        taskLabels{i} = taskLabels;
        seedSigns{i}   = seedSignLabels(CAPix{i});

         
         % Temporal normalization, as done by Karahanoglu et al., 2015
        if strcmp(dist, 'cosine')
            sd1 = std(cluster{i}, 0, 2);
            if (~ismember(0, sd1(voxels)) && ~ismember(Inf, sd1(voxels)))
                CAPmean(voxels) =  bsxfun(@rdivide,  CAPmean(voxels), sd1(voxels));
                CAPmean{i}      = CAPmean;
            else
                warning(['Voxels were not normalised by variance for cap=' num2str(i) ', K=' num2str(k)]);
            end
        end
        
        % Spatial nomalization (voxelwise, for visualisation purposes only)
        sd2 = std(CAPmean(voxels));
        if (~ismember(0, sd2) && ~ismember(Inf, sd2))
            CAPmean(voxels)       = (CAPmean(voxels) - mean(CAPmean(voxels)))/sd2;
            CAPmean{i}     = CAPmean;
        end
        
        CAP_vol{i}     = reshape(CAPmean{i} , 53, 63, 46); % save it back in the original 3D shape
   
        
        % Creates a new NIFTI for the CAP
        tmp_cap = make_nii(CAP_vol{i},voxel_size,round(-voxel_shift./voxel_size));
        tmp_cap.hdr.dime.datatype=64;
        tmp_cap.hdr.dime.bitpix=64;
        
        % Saves the cap nifti
        capNiftiFile = fullfile(ppicapsDir,['ppicaps_2Level_cap' num2str(i) '_' date '.nii']);
        save_nii(tmp_cap,capNiftiFile);
      
end

end


% From here on, statistical analysis is performed using the same code as in
% that used in ppicaps_toyExample.m for the statistical analysis part.




% ------------------------------------------------------
% HELPER FUNCTIONS
% ------------------------------------------------------

% Finds indeces of voxels in or out of Grey Matter (GM) 
function [inGM, nonGM] = findGMvoxels

mask      = load([pwd '/dependencies/MNIBrainMask.mat']); % Load MNIBrainMask
nonGM     = find(mask.MNIBrainMask(:) < 1); voxels  = find(mask.MNIBrainMask(:) == 1);

end

function [voxel_size, voxel_shift] = prepareToSaveNifti

VOI = load([pwd '/dependencies/VOI_PCC_1.mat']);
voxel_size = diag(VOI.xY.spec.mat);
voxel_size = voxel_size(1:end-1)';

voxel_shift = VOI.xY.spec.mat(:,4);
voxel_shift = voxel_shift(1:end-1)';

end