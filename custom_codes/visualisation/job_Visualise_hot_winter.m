% JOB_Visualise: This function plots results from a subject from
% the Movie Watching fMRI dataset using the SLOVER function.
%____________________________________________________________________________
% Copyright (C) 2016 MIP:Lab

% Lorena Freitas
% $Id: job_Visualise.m 11 2016-12-10 20:10:27 F Lorena $
%
% Modified by Serafeim Loukas, July 19, 2021

function img = job_Visualise_hot_winter(capNiftiFile, frame)

if isempty(frame)
    frame = 'axial'; % default view
end

%% set paths
codeBasePath = '/Users/loukas/Desktop/PPI_CAPs/scripts/visualisation/';
%structFile = '/Users/loukas/Desktop/PPI_CAPs/template/template_T2_fused40w_n20.nii';
%structFile='/Users/loukas/Desktop/PPI_CAPs/subjects/PTMT1_49RiGo_Rechercheneonat_3151/structural/wT2.nii';
structFile='/Users/loukas/Desktop/PPI_CAPs/template/spm152_wrapped2neonatal.nii';

%% choose case
switch frame
    case 'axial'
        obj = load(fullfile(codeBasePath,'three_layers_object.mat'));
        obj = obj.o;
    case 'sagittal'
        obj = load(fullfile(codeBasePath,'three_layers_object.mat'));
        obj = obj.o; obj.transform = 'sagittal';
    case 'coronal'
        obj = load(fullfile(codeBasePath,'three_layers_object.mat'));
        obj = obj.o; obj.transform = 'coronal';
end

set(gcf,'Position',[484 77 560 751]);
obj.figure=gcf;

% this object has three layers, the structural one and two overlays (one
% for positive and one for negative values)
obj.img(1).vol=spm_vol(structFile); % file names
obj.img(2).vol=spm_vol(capNiftiFile);
obj.img(3).vol=spm_vol(capNiftiFile);

% get min and max values of layers
[mx1 mn1] = slover('volmaxmin', obj.img(1).vol);
[mx mn] = slover('volmaxmin', obj.img(2).vol);

% thresholds for layer 2 and 3
lowerthreshold = 0.05; 
%upperthreshold = mx;
upperthreshold = max(abs(mx) , abs(mn)) / 2;  % custom max value for nice visualizations

% define color map range
obj.img(1).range = [35 85]; % [35 85] for spm152_wrapped2neonatal.nii as structural / [0 mx1] for others
obj.img(2).range=[lowerthreshold, upperthreshold];    % make symmetric
obj.img(3).range=[-lowerthreshold, -upperthreshold];  % make symmetric

% define color maps
obj.img(1).cmap = obj.img(1).cmap * 1.1; % make struct. darker (<1) or brighter (>1)
obj.img(2).cmap='hot';
obj.img(3).cmap='winter'; 

% transparency % type
obj.img(1).type = 'truecolour';
obj.img(2).type = 'split';
obj.img(3).type = 'split';
obj.img(2).prop = 1;
obj.img(3).prop = 1;

% define which slices to plot
if strcmp(frame,'coronal')
   obj.slices = -50:3:15;
elseif strcmp(frame,'sagittal')
    obj.slices = -35:3:30;
else
   obj.slices=-22:3:30;
end

% plot updated figure
paint(obj);

% assign output variable
img = gcf;

end



% % this object has two layers, the structural one and one overlay for the
% % contrast
% structural = dir(['/Volumes/EPFL_Lorena/Autistes/Preproc_ArtRepair' '/MNI*.nii']);
% obj.img(1).vol=spm_vol(['/Users/lorenafreitas/Software and Tools/spm12/canonical/single_subj_T1.nii']); % structural
%
% % for the second layer, I load an arbitrary file to fill the struct fields
% % that I don't care about, then I update just the "fname" field
% %files=spm_select('FPListRec',deblank(b.dir_fonc(1,:)),['^s6w3.*\.img$']);
% %ResultVOL = spm_vol('/Volumes/EPFL_Lorena/Autistes/Preproc_ArtRepair/PPI_05-Dec-2016_Level2/PPI_STSx(Fun-Science)_05-Dec-2016_Level2/All10Subjects/spmT_0001.nii');
% %ResultVOL = spm_vol(['/Volumes/EPFL_Lorena/Autistes/Preproc_ArtRepair/PPI-CAPs_2Level_sMap_', date, '.nii']);
% %ResultVOL = spm_vol(['/Volumes/EPFL_Lorena/Autistes/Preproc_ArtRepair/sMap_T2_norm_mask.nii']);
% ResultVOL = spm_vol(capNiftiFile);
%
% %ResultVOL.fname = [b.dataDir b.curSubj '_deconvolvedData'];
% %ResultVOL.fname = '/Volumes/EPFL_Lorena/Autistes/Preproc_ArtRepair/PPI_05-Dec-2016_Level2/PPI_STSx(Fun-Science)_05-Dec-2016_Level2/All10Subjects/spmT_0001.nii';
% % //TODO: CREATE .nii FOR DECONVOLVED DATA IN ORDER TO BE ABLE TO LAYER IT IN SLOVER
% %ResultVOL = spm_vol('/Users/lorenafreitas/Software and Tools/ndhist/zerosinPPI.nii');
%
% obj.img(2).vol=ResultVOL;%spm_vol([newestGLM '/spmT_0002.nii']); % Contrast Sound from GLM (2 = Sounds / 7 = Mother-Stranger)
%
% % get min and max values of layers
% [mx mn] = slover('volmaxmin', obj.img(2).vol);
% a = max(abs(mx), abs(mn));
% %mx = a;
% %mn = -a;
% %mn = -mx;
% mn =1;
% %mn = -1; mx =1;
% % define color map range
% obj.img(2).range=[mn  mx];
%
% % define color maps
% obj.img(2).cmap='hot';
% %obj.img(1).cmap=obj.img(1).cmap*1.5; % making struct brighter
%
% % define which slices to plot
% obj.slices=[24];
% %obj.slices=-40:2:68;
% %obj.slices = 30;
%
% % plot updated figure
% obj = paint(obj);
%
% img = gcf;
% % save figure in a folder where it's easy to analyse
% %saveas(gcf,sprintf('%s/myPLS_%s_BS%d_NORM%d%d_brain',CONST_OUTPUT_PATH,CONST_CONDITION,iter_lv,CONST_NORM_IMAGING,CONST_NORM_BEHAV),'epsc');
%
% end