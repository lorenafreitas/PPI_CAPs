function step_3_main_PPI_CAPs_sl
% -------------------------------------------
%    Author: Serafeim Loukas, July 13, 2021
%    For PPI-CAPs project with Joana (& Lorena)
%
%    This is the main script for the PPI-CAPs analysis assuming that 
%    the 2 previous steps have been sucessfully executed 
%    (seed creation & subselection of the frames).
%
%    Updated & modified by: Serafeim Loukas on Oct 27, 2021
% -------------------------------------------
% main_PPI_CAPs_sl - function to execute the main PPI-CAPs workflow

clc;
clear all;

%% Set seed - reproduce results
% Note: in some rare cases, PPI-CAPs might look slightly different when re-running the code
% k-means with 100 replicates is stable but it can also fall in local minima...
rng('default')
rng(6)

%% General paths setup
% set main path (universal path)
% -------------------------------
run_consensus = 0; % run consensus clustering 1=yes, 0=no (if already ran)
run_clustering = 0; % 1=run clustering of frames, 0=load saved results (if already ran)
run_stats = 1; % perform stat analysis (always 1 in principal)

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

%% T2 template (this step is not needed)
% -------------------------------
%T2_template_file = [ppicapsDir, 'template/', 'template_T2_fused40w_n20.nii'];

%% Create cell with the subject folder names
% -------------------------------
n_subjects= size(S,1);
subject_name= {};
for k = 1 : n_subjects
    subject_name{k} = subFolders(k).name;
end

%% Load saved data for all subjects (this matrix was produced using the "main_threshold_frames_sl.m" script)
% change this depending the case e.g. this is the 5% retained frames case
% -------------------------------
% load as global variables to avoid loading them twice (for consensus clustering) since they are huge
global PPIframesALL; global subjectLabelALL;
file_to_load = 'PPIallInfo_5.mat';
disp(['Loading the data: ' file_to_load])
load([ppicapsDir file_to_load]);
disp('Data loaded')

%% Find voxels in and out of grey matter (global grey matter mask for all subjects)
% -------------------------------
[voxels, nonVoxels] = ppicaps_util_getMaskedVoxels(ppicapsDir);

%% Print some info
disp(['The number of GM voxels is: ' num2str(size(voxels,1))])
disp(['The number of group-level frames is: ' num2str(size(PPIframesALL,2))])

%% Find matrix rotations for saving nifti files
% -------------------------------
[voxel_size, voxel_shift] = prepareToSaveNifti(ppicapsDir);

%% Run consensus clustering to determine best K
% -------------------------------
disp(['Starting the consensus clustering for data matrix: ' file_to_load])
% get case name
mytmp = strsplit(file_to_load, '_');mytmp2 = mytmp{2};mytmp3 = strsplit(mytmp2, '.'); % get case name

% run the consensus Clustering
if run_consensus
    main_consensusClustering(ppicapsDir, file_to_load, mytmp3{1})
    close all % close open figures (figs have been already saved)
end
% NOTE: After getting the plots, choose the best K and define it below for the main clustering.

%% Cluster timepoints, using frames as feature vectors
% Choose number of clusters (PPI-CAPs) -- initial guess (or based on the literature)
% -------------------------------
k = 4; % initial guess or based on the literature -- should be defined using consensus clustering
%k = input("Define best K: ");

%% Distance used for kmeans 
% -------------------------------
dist = 'cosine';
replicates = 100; % k-means # of internal replicates to run

%% Run kmeans++ on all selected suprathresholded frames
% PPIframesALL : group-level whole brain voxel-wise frames
% -------------------------------
if run_clustering
    disp(['Starting kmeans clustering with k=' num2str(k)]);
    [PPI_CAPs, PPI_CAPs_Ix]  = kmeanspp(PPIframesALL(voxels,:), k, dist, replicates); % "PPIframesALL(voxels,:)": use only GM voxels for clustering
    disp('Kmeans clustering is done !')
    %save the clusters so that it is not needed to rerun for new contrasts
    save([save_caps 'PPI_CAPs5'], 'PPI_CAPs'); save([save_caps 'PPI_CAPs_Ix5'], 'PPI_CAPs_Ix');
else % load saved results
    disp('Not running clustering - loading already saved results');
    load([save_caps 'PPI_CAPs5']); load([save_caps 'PPI_CAPs_Ix5']);
end

%% main loop, for each estimated CAP/cluster
% -------------------------------
disp('Working on the estimated CAPs')
for i = 1:k
        CAPix{i}              = find(PPI_CAPs == i); % #frames x 1; Indices of frames for this CAP
        CAPl{i}               = PPI_CAPs == i; % boolean mask of frames for this CAP
        scf{i}                = CAPl{i}.* PPI_CAPs_Ix; % signed cluster's frames
        cluster{i}            = PPIframesALL(:,CAPix{i}); % #voxels x #frames, get only the frames that belong to this clsuter
        CAPmean{i}            = (sum(PPIframesALL(:,scf{i}==1),2) - sum(PPIframesALL(:,scf{i}==2),2)) / sum(CAPl{i}); % whole-brain centroid (average)
        assert(size(CAPix{i},2)==sum(CAPl{i})); % sanity check 
        CAPmean{i}(nonVoxels,:)  = 0; % force zero-out non-brain voxels
        
        polarityLabelsTmp= scf{i}; % polarities of all group-level frames
        polarityLabels{i} = polarityLabelsTmp(CAPix{i}); % polarities only for the frames of the current cluster
        assert(isequal(polarityLabelsTmp(CAPix{i}) , polarityLabelsTmp(CAPl{i}))); % sanity check 
        taskLabels{i} = task4framesALL(CAPix{i}); % task encoding variable for the frames that belong to the current cluster
        seedSigns{i}   = seedSignALL(CAPix{i}); % seed sign encoding variable for the frames that belong to the current cluster
        
        % Temporal normalization of the centroid, as done by Karahanoglu et al., 2015
        if strcmp(dist, 'cosine')
            sd1 = std(cluster{i}, 0, 2); % vector: std of each voxel estimated across frames (rows)
            CAPmean{i}(sd1(:)>0) =  bsxfun(@rdivide,  CAPmean{i}(sd1(:)>0), sd1(sd1(:)>0)); % normalize centroid by std
        end
        
        % Spatial nomalization - within-brain z-scoring (voxelwise, for visualisation purposes only)
        sd2 = std(CAPmean{i}(CAPmean{i}~=0)); % this is a scalar i.e. the std of the centroid/cluster
        CAPmean{i}(CAPmean{i}~=0)       = (CAPmean{i}(CAPmean{i}~=0) - mean(CAPmean{i}(CAPmean{i}~=0))) / sd2; % norm only within-brain
        
        % Reshape back to 3D
        CAP_vol{i}     = reshape(CAPmean{i} , 79, 95, 68); % save it back in the original 3D shape -> mudado numeros para dimensoes dos volume functional normalisado rtemplate
        
        % Creates a new NIFTI for the CAP
        mytmp = strsplit(file_to_load, '_');mytmp2 = mytmp{2};mytmp3 = strsplit(mytmp2, '.');
        %capNiftiFile = fullfile(save_caps, ['PPI-CAPs_2Level_K' num2str(k) '_cap' num2str(i) '_' date '.nii']); 
        capNiftiFile = fullfile(save_caps, ['PPI-CAPs_2Level_K' num2str(k) '_cap' num2str(i) '_thres_' mytmp3{1} 'perc' '.nii']); 
        
        % save using save_nii - this does not result in realigned maps.......
        %tmp_cap = make_nii(CAP_vol{i}, voxel_size, round(-voxel_shift./voxel_size));
        %tmp_cap.hdr.dime.datatype=64;
        %tmp_cap.hdr.dime.bitpix=64;
        % Saves the cap nifti
        %save_nii(tmp_cap, capNiftiFile);
        
        % save using spm because the other way does not work as expected
        header = spm_vol([ppicapsDir 'subjects/' 'TemplateGMmask.nii']);
        header.fname = capNiftiFile;
        header.descrip = '';header.pinfo = [0;0;352];% force re-scaling
        spm_write_vol(header, CAP_vol{i});
end

%% -------------------------------------------------------------------------
% STATISTICAL ANALYSIS
% -------------------------------------------------------------------------
if ~run_stats; disp('Will not run statistical analysis');else; disp('Starting statistical analysis'); end

if run_stats
    % Confusion matrices and permutation testing to check whether each 
    % PPI-CAP's frames correlate with any of the effects we are trying to analyse.
    seedSigns = seedSignALL;
    taskSigns = task4framesALL;
    grouplabels = groupLabelALL;
    
    %% Build contrast of interest e.g. music vs singing
    % 1: flute, 2: mother, 3: stranger, 4: noise, 5: silence
    % unique(task4framesALL) % prints 0,1,2,3,4,5
    
    % Contrast: music VS singing (combined)
    task_mask = ~(taskSigns==1 | taskSigns==2 |  taskSigns==3);
    taskSigns(task_mask)=0; % zero-out all except 1,2,3 that interests us
    taskSigns(taskSigns==1) = 1;  % as in the original publication: music VS 
    taskSigns(taskSigns==2) = -1; % as in the original publication: stranger + mother singing
    taskSigns(taskSigns==3) = -1; % as in the original publication: stranger will be combined with mother singing
    contrast_name = 'musicVSsinging';

%     % Contrast: music VS silence
%     task_mask = ~(taskSigns==1 | taskSigns==5);
%     taskSigns(task_mask)=0; % zero-out all except 1,5 that interests us
%     taskSigns(taskSigns==1) = 1;  % as in the original publication: music VS 
%     taskSigns(taskSigns==5) = -1; % as in the original publication: silence
%     contrast_name = 'musicVSsilence';

    
    %% Build group contrast of interest e.g. PTC vs PTM
    % 1:PTC, 2:PTM, 3:FT
%     group_mask = ~(grouplabels==1 | grouplabels==2);
%     grouplabels(group_mask)=0;
%     grouplabels(grouplabels==1) = 1; % PTC
%     grouplabels(grouplabels==2) = -1; % PTM
    
    %% Permutation & plotting params
    maxPerm = 2000; % Maximum number of random permutations
    f=figure; f.Position(3:4) = [1000 850];
    colors = cbrewer('seq', 'YlOrRd', 50, 'pchip');

    %% Calculate confusion matrices for each PPI-CAP.
    for i = 1:k
        thisCapSeed = seedSigns(PPI_CAPs==i);
        thisCapTask = taskSigns(PPI_CAPs==i);
        thisCapPPI = thisCapSeed .* thisCapTask; % element-by-element multiplication
        thisCapClusteringFlip = PPI_CAPs_Ix(PPI_CAPs==i);
%         thisCapGroupLabels = grouplabels(PPI_CAPs==i);

        % Seed
        subplot(k, 3, 3*i-2);
        [seedConfusionMatrix{i}, seed_fkix{i}] =  calc_confusionMatrix(thisCapSeed, thisCapClusteringFlip);
        imagesc(seedConfusionMatrix{i}); colormap(colors);colorbar;

        % Task
        subplot(k, 3, 3*i-1);
        [taskConfusionMatrix{i},  task_fkix{i}] = calc_confusionMatrix(thisCapTask, thisCapClusteringFlip);
        imagesc(taskConfusionMatrix{i});colormap(colors);colorbar;

        % PPI
        subplot(k, 3, 3*i);
        [ppiConfusionMatrix{i}, ppi_fkix{i}]  = calc_confusionMatrix(thisCapPPI, thisCapClusteringFlip);
        imagesc(ppiConfusionMatrix{i});colormap(colors);colorbar;
        
%         % Group
%         subplot(k, 3, 3*i);
%         [ppiConfusionMatrix{i}, ppi_fkix{i}]  = calc_confusionMatrix(thisCapGroupLabels, thisCapClusteringFlip);
%         imagesc(ppiConfusionMatrix{i});colormap(colors);colorbar;        

        % Random Permutations
        for nperm = 1: maxPerm
            % --- seed
            randThisCapSeed = thisCapSeed(randperm(length(thisCapSeed)));
            [~, seed_fkix_rnd(i,nperm)]= calc_confusionMatrix(randThisCapSeed, thisCapClusteringFlip);   
            % --- task
            randThisCapTask = thisCapTask(randperm(length(thisCapTask)));
            [~, task_fkix_rnd(i,nperm)]= calc_confusionMatrix(randThisCapTask, thisCapClusteringFlip);   
            % --- ppi
            randThisCapPPI = thisCapPPI(randperm(length(thisCapPPI)));
            [~, ppi_fkix_rnd(i,nperm)]= calc_confusionMatrix(randThisCapPPI, thisCapClusteringFlip);
        end
    end

    %% Now let us plot the distributions of fk-indeces obtained through random 
    % permutations of the (seed; task; and ppi) labels from the frames that
    % make up each PPI-CAP. 
    % Then, we plot a line for the real data-driven fk-index and see where it stands!
    
    f=figure; f.Position(3:4) = [1000 850];
    for i = 1:k     
        % Seed
        ax2 = subplot(k, 3, 3*i-2);
        h2 = histogram(ax2,seed_fkix_rnd(i,:));title(ax2, 'Null Seed Distribution');
        h2.BinMethod = 'Sturges'; h2.Normalization = 'countdensity'; h2.FaceColor = colors(i,:); h2.NumBins = 10; set(gca,'FontSize',11);
        line(ax2, [seed_fkix{i}, seed_fkix{i}], get(ax2, 'ylim'), 'LineWidth', 2.5, 'Color', 'r');
        pSeed = min(((length(find(seed_fkix_rnd(i,:)<seed_fkix{i}))+1)/(maxPerm+1)),((length(find(seed_fkix_rnd(i,:)>seed_fkix{i}))+1)/(maxPerm+1)));
        text(ax2, 0.7,0.9,['p = ' num2str(pSeed)],'Units','normalized','FontSize',13);
        tmp_ = gca;
        tmp_.TitleFontSizeMultiplier = 1.3;
        
        % Task
        ax4 = subplot(k, 3, 3*i-1);
        h2 = histogram(ax4,task_fkix_rnd(i,:));title(ax4, 'Null Task Distribution');
        h2.BinMethod = 'Sturges'; h2.Normalization = 'countdensity'; h2.FaceColor = colors(i,:);h2.NumBins = 10; set(gca,'FontSize',11);
        line(ax4, [task_fkix{i}, task_fkix{i}], get(ax4, 'ylim'), 'LineWidth', 2.5, 'Color', 'r');
        pTask = min(((length(find(task_fkix_rnd(i,:)< task_fkix{i}))+1)/(maxPerm+1)),((length(find(task_fkix_rnd(i,:)>task_fkix{i}))+1)/(maxPerm+1)));
        text(ax4, 0.7,0.9,['p = ' num2str(pTask)],'Units','normalized','FontSize',13);
        tmp_ = gca;
        tmp_.TitleFontSizeMultiplier = 1.3;
        
        % PPI
        ax3 = subplot(k, 3, 3*i);
        h2 = histogram(ax3, ppi_fkix_rnd(i,:));title(ax3, 'Null PPI Distribution');
        h2.BinMethod = 'Sturges'; h2.Normalization = 'countdensity'; h2.FaceColor = colors(i,:);h2.NumBins = 10; set(gca,'FontSize',11);
        line(ax3, [ppi_fkix{i}, ppi_fkix{i}], get(ax3, 'ylim'), 'LineWidth', 2.5, 'Color', 'r');
        pPPI = min(((length(find(ppi_fkix_rnd(i,:)< ppi_fkix{i}))+1)/(maxPerm+1)),((length(find(ppi_fkix_rnd(i,:)>ppi_fkix{i}))+1)/(maxPerm+1)));
        text(ax3, 0.7,0.9,['p = ' num2str(pPPI)],'Units','normalized','FontSize',13);
        tmp_ = gca;
        tmp_.TitleFontSizeMultiplier = 1.3;
        
    end
    
    %% Save results as .png figures
    h = get(0,'children');
    for i=1:length(h)
        print(h(i), fullfile(save_caps,['stat_figure_' num2str(i) '_thres_' mytmp3{1} 'perc_' contrast_name]), '-dpng', '-r300','-painters');
    end
    
end
disp('All done !!!')
end

