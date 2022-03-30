% -------------------------------------------------------
%
%    thresholdFrames - function to threshold frames
%
%    Created by:         Lorena Freitas
%    Last checked:       28.09.2019
%
%    Updated & modified by: Serafeim Loukas on Sep 10, 2021
%
%    BUG FIX: z-scoring across time and not frame-wise (line 36)
%    BUG FIX 2: z-scoring each voxel across time (line 109)
%    
%    Updated & modified by: Serafeim Loukas on Sep 7, 2021
%
%    Point VOI to universal location instead of pointing to each subject's subfolder (faster, lines 49-51) 
% ------------------------------------------------------

function [PPIframes, thresholdedIx, task4frames, seedSign] = thresholdFrames(b, seed, percent)

% ----------------------
% Variable set up
% ----------------------

% SPM's seed name
VOI_Name    = ['VOI_' seed '_1'];

% Folder where analysis files were stored
%b           = initialize_vars(thisSubject);
analyze_dir = [b.dataDir];

% Load subject's data (deconvData) that are in the functional space
load(spm_select('FPListRec',[b.dataDir],'_deconvolvedData.mat$'));

% Total number of frames
nFrames     = size(deconvData,4); % N of "good" frames because deconvolution was done only using the good volumes
disp(['good vols:' num2str(nFrames)])

% Reshape matrix into 2D to facilitate computations & z-score each voxel across time (frames)
% FIXED: z-score each voxel across time
%subjectData2D      = zscore(reshape(deconvData, [], nFrames))';
subjectData2D       = zscore(reshape(deconvData, [], nFrames), 0, 2)'; % FIXED: z-score each voxel across time/frames

% --------------------------------------------
% Create seed timecourse (VOI creation was done in step_1 script
% --------------------------------------------

% Obtain seed, seed xyz locations and extract mean signal from seed
%VOI         = spm_select('FPListRec',[analyze_dir 'results_mask_new' filesep 'Voice_fMRI_2019_RPARAM24'],['^' VOI_Name '.mat$']); VOI = load(deblank(VOI(1,:)));

regtmp = regexp(b.dataDir,filesep,'split');
seed_path = [filesep fullfile(regtmp{1:6}) filesep];
VOI         = spm_select('FPListRec',seed_path,['^' VOI_Name '.mat$']); VOI = load(deblank(VOI(1,:)));
VOIxyzMNI   = VOI.xY.XYZmm; % Dimensions: [3 x #voxels double]
disp(['the # of voxels in VOI (func space) is:' num2str(size(VOI.xY.XYZmm,2))])
VOIxyzMat   = round(VOI.xY.spec.mat \ [VOIxyzMNI; ones(1, length(VOI.xY.XYZmm))]);
VOIxyzMat   = VOIxyzMat(1:3,:); % XYZ indeces of seed's voxels in "MATLAB matrix" space
clear VOI_Name VOI VOIxyzMNI

% Average the seed signal from all voxels
meanSeedSignal = averageSignal(deconvData, VOIxyzMat, nFrames); % using only good volumes

% ----------------------
% Create task timecourse
% ----------------------

% Repetition time
TR = 0.7;

% Onset times for task FLUTE  (fun scences), in seconds
onset_FLUTE = [8.307 65.184 89.531 130.072 186.894 227.475 251.796 341.022 365.384 397.846];

% Onset times for task MOTHER (fun scences), in seconds
onset_MOTHER = [16.432 40.818 73.302 154.430 178.771 211.247 259.900 284.236 349.153 373.503];

% Onset times for task STRANGER (fun scences), in seconds
onset_STRANGER = [24.563 48.948 105.736 146.306 195.030 219.363 276.122 308.557 332.897 381.625];

% Onset times for task NOISE (fun scences), in seconds
onset_NOISE = [32.693 81.408 113.841 138.175 170.648 203.142 235.579 300.446 324.766 357.277];

% Onset times for task SILENCE (fun scences), in seconds
onset_SILENCE = [0.180 57.061 97.629 121.947 162.535 243.687 268.009 292.340 316.660 389.734 405.948];

% All onsets
onsets = {onset_FLUTE, onset_MOTHER, onset_STRANGER, onset_NOISE, onset_SILENCE};

% Duration of blocks. 
dur = 8;

SPM = spm_select('FPListRec',analyze_dir,'^SPM.mat$'); load(SPM);
BAD = spm_select('FPListRec',[analyze_dir 'functional/'],'^motion_0.5.mat$'); load(BAD);

% Get task labels for all frames
[taskLabels]    = createTaskRegressor(onsets, dur, TR, sess);

% Find index of frames to threshold
[thresholdedIx, ~, seedSign] = thresholdSeed(meanSeedSignal, nFrames, percent);

% Find suprathreshold frames to keep (selected frames for the current subject)
PPIframes     = subjectData2D(thresholdedIx,:);
disp(['selected vols:' num2str(size(PPIframes,1))])

% Find task labels for suprathreshold frames (i.e., which task was performed when the selected frames were acquired?)
task4frames   = taskLabels(thresholdedIx);

end


% ----------------
% HELPER FUNCTIONS
% ----------------

function meanSignal = averageSignal(X, VOIxyzMat,nFrames)

% Initialise variables
nVoxels = size(VOIxyzMat,2);
signal  = zeros(nVoxels, nFrames);
%X       = zscore(X);
X       = zscore(X, 0, 4); % z-score each voxel across time ! e.g. mean(X(20,20,20,:)), std(X(20,20,20,:)), plot(squeeze(X(20,20,20,:)))

% Get signal from voxels
for i=1:nVoxels
    x = VOIxyzMat(1,i); y = VOIxyzMat(2,i); z = VOIxyzMat(3,i);
    signal(i,:) = squeeze(X(x,y,z,:))';
end

% Average signal from voxels
meanSignal = mean(signal); % meanSeedSignal = [1 x # nscan]

end


function [taskLabels] = createTaskRegressor(onsets, dur, TR, sess)
% change this as well
nscans = 590;
TMP_timecourse = zeros(1,nscans);
nSessions = size(sess,1);

for thisSession = 1:nSessions
    TMP_timecourse(sess(thisSession,1):sess(thisSession,2)) = 1; % data were split in some case to make sessions due to motion
end

thisDuration = round(dur/TR); % # scans
taskLabels    = zeros(1,nscans); 

for thisTask = 1:length(onsets)
    scanTimeOnsets = ceil((onsets{thisTask}/TR)); % in scans
    for thisOnset = 1:length(scanTimeOnsets)
        taskLabels(scanTimeOnsets(thisOnset):(scanTimeOnsets(thisOnset)+thisDuration-1)) = thisTask; % 1: flute, 2: mother, 3: stranger, 4: noise, 5: silence
    end
end

taskLabels(TMP_timecourse==0)=[]; % motion corrupted timepoints (scans) outside the sess limits are set to 0

end


function [thresholdedIx, threshold, seedSign] = thresholdSeed(meanSeedSignal, nFrames, percent) 

nFrames2keep                            = floor(nFrames*percent / 100);
seedSign                                = zeros(1, nFrames2keep); 
[sortedSeedS,sortingIxS]                = sort(meanSeedSignal,'descend');
PosFrames2keep                          = round(nFrames2keep/2); % half will be positive seed moments
NegFrames2keep                          = nFrames2keep - PosFrames2keep; % the other half will be negative seed moments
thresholdedIx                           = sortingIxS(1:PosFrames2keep);
threshold                               = min(meanSeedSignal(thresholdedIx));
seedSign(1:PosFrames2keep)              = 1;
seedSign(PosFrames2keep+1:end)          = -1;

% Now take also moments where the seed is deactive !
thresholdedIxNeg            = sortingIxS((end - NegFrames2keep + 1):end); % index of most extreme negatives
[thresholdedIx, sortedIxIx] = sort([thresholdedIx thresholdedIxNeg]);

seedSign = seedSign(sortedIxIx);

end

