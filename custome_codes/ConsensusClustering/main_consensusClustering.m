function main_consensusClustering(dataDir, file_to_load, case_name)
% -------------------------------------------
% This script uses ConsensusClustering to find the best number of Ks
%    Lorena Freitas
%    $Id: ppicaps_util_consensusClustering.m 11 2018-11-27 17:26:24F Lorena $
%
%    Updated & modified by Serafeim Loukas on July 19, 2021
%     - Made "outDir" case-specific to avoid conflict when running the same script in parallel (e.g. for other
%    thersholds)
%    
% -------------------------------------------

%% Load data for classification
% -------------------------------
% % base data directory
% disp(['Loading the data: ' file_to_load])
% load([dataDir file_to_load]);
% disp('Data loaded')

% instead of re-loading use global variables
global PPIframesALL; global subjectLabelALL;

% The group-level data
X = PPIframesALL'; 
subject_labels = subjectLabelALL';
%clear PPIframesALL subjectLabelALL %taskInfo timeInfo

%% Initialise parameters
% -------------------------------
highlyActiveVoxelsOnly = 0;
param.K = 3:8;
param.cons_n_folds = 10; % number of subsamples for each k
param.Subsample_fraction = 0.8;
param.Subsample_type = 'subjects';
param.DistType = 'cosine';
%outDir = [dataDir 'results/consensusClustering/'];
outDir = [dataDir 'results/' 'consensusClustering_' case_name '/'];
param.KmeansMethod = 'kmeanspp'; % kmeansmatlab or kmeanspp

%% Mask non-GM voxels
% -------------------------------
[voxels, nonVoxels] = ppicaps_util_getMaskedVoxels(dataDir);
consensus_input = X(:, voxels); % frames (all subjects) x only GM voxels
disp(['The number of GM voxels are: ' num2str(size(voxels,1))])

%% Use only highly active voxels
% -------------------------------
if highlyActiveVoxelsOnly
    X = ppicaps_util_highSNRmask(X,6);
end

%% Run consensus clustering based only on GM voxels (see: "consensus_input" variable)
% -------------------------------
result = ConsensusClustering(consensus_input, subject_labels, param, outDir, case_name);

%% Plot consensus clustering results (better run this locally)
% -------------------------------
consensusClusteringPlotResults(result, outDir, case_name);
end