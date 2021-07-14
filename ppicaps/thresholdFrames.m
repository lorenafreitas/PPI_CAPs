% -------------------------------------------------------
%
%    thresholdFrames - function to threshold frames
%
%    Created by:         Lorena Freitas
%    Last checked:       28.09.2019
%
%
% ------------------------------------------------------

function [PPIframes, thresholdedIx, task4frames] = thresholdFrames(thisSubject, seed, percent)

% ----------------------
% Variables' set up
% ----------------------

% SPM's seed name
VOI_Name    = ['VOI_' seed '_1'];

% Folder where analysis files were stored
b           = initialize_vars(thisSubject);
analyze_dir = [b.dataDir 'analysis/'];

% Total number of frames
nFrames     = 179;

% Load subject's data (deconvData)
load([b.dataDir b.curSubj '_deconvolvedData']);

% Reshape matrix into 2D to facilitate computations
subjectData2D      = zscore(reshape(deconvData, [], nFrames), 0, 2)'; % z-score each voxel across frames/time


% ----------------------
% Create seed timecourse
% ----------------------

% Obtain seed
VOI         = spm_select('FPListRec',analyze_dir,[VOI_Name '.mat$']); VOI = load(VOI);
VOIxyzMNI   = VOI.xY.XYZmm; % Dimensions: [3 x #voxels double]
VOIxyzMat   = round(VOI.xY.spec.mat \ [VOIxyzMNI; ones(1, length(VOI.xY.XYZmm))]);
VOIxyzMat   = VOIxyzMat(1:3,:); % XYZ indeces of seed's voxels in "MATLAB matrix" space
clear VOI_Name VOI VOIxyzMNI

% Average the seed signal from all voxels
meanSeedSignal = averageSignal(deconvData, VOIxyzMat);

% ----------------------
% Create task timecourse
% ----------------------

% Repetition time
TR_MW     = 2;

% Onset times for task A (fun scences), in seconds
onset_Fun = [5 22 89 182 288];

% Duration of blockes. First row: dur_fun; Second row: dur_sci
dur_MW    = [13 8 45 52 23; 4 59 48 54 42];

% Get task labels for all frames
[~, taskLabels]    = createTaskRegressor(onset_Fun, dur_MW, TR_MW, nScans);

% Find index of frames to threshold
[thresholdedIx, ~] = thresholdSeed(meanSeedSignal, nFrames, percent);

% Find suprathreshold frames to keep
PPIframes     = subjectData2D(thresholdedIx,:);

% Find task labels for suprathreshold frames
task4frames   = taskLabels(thresholdedIx);

end





% ----------------
% HELPER FUNCTIONS
% ----------------


function meanSignal = averageSignal(X, VOIxyzMat)

% Initialise variables
nVoxels = size(VOIxyzMat,2);
signal  = zeros(nVoxels, 179);
X       = zscore(X);

% Get signal from voxels
for i=1:nVoxels
    x = VOIxyzMat(1,i); y = VOIxyzMat(2,i); z = VOIxyzMat(3,i);
    signal(i,:) = squeeze(X(x,y,z,:))';
end

% Average signal from voxels
meanSignal = mean(signal); % meanSeedSignal = [1 x # nscan]

end


function [taskRegressor, taskLabels] = createTaskRegressor(onset_Fun, dur_MW, TR_MW, nScans)

scanTimeFunOnsets = ceil((onset_Fun/TR_MW)); % round up (?)
taskRegressor1    = zeros(1,nscans);

thisOnset    = scanTimeFunOnsets(i);
thisDuration = round(dur_MW(1,i)/TR_MW);
taskRegressor1(thisOnset:(thisOnset+thisDuration)) = 1;

% Flip zeros and ones
taskRegressor2 = 1 - taskRegressor1; 

% After the 178th scan it's just resting state
taskRegressor2([1,2,175:end]) = 0; 
taskRegressor = taskRegressor1;

taskLabels = taskRegressor;
taskLabels(taskRegressor2==1) = 2;
 % create task regressor with 1 / -1 values
taskRegressor(taskRegressor2==1) = -1;

end



function [thresholdedIx, threshold] = thresholdSeed(meanSeedSignal, nFrames, percent)

nFrames2keep              = floor(nFrames*percent / 100);
[sortedSeedS,sortingIxS]  = sort(meanSeedSignal,'descend');
thresholdedIx             = sortingIxS(1:ceil(nFrames2keep/2));
threshold                 = min(meanSeedSignal(thresholdedIx));

% Now take also moments where the seed is deactive
thresholdedIx            = sort([thresholdedIx sortingIxS(end-ceil(nFrames2keep/2):end)])

end