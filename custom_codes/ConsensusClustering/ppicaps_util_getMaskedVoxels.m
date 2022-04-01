function [voxels, nonVoxels] = ppicaps_util_getMaskedVoxels(path)

% Exclude non-brain voxels
% the GM mask was resliced to TemplateGMmask.mat (resliced in the functional space)
TemplateGMmask = load([path 'subjects/' 'TemplateGMmask.mat']); %mudar mascara GM da template 40w -> Vout (TemplateGMmask.mat)

% find voxels inside / outside of the brain
nonVoxels = find(TemplateGMmask.maskMappedToFmriSpace(:) < 1); 
voxels  = find(TemplateGMmask.maskMappedToFmriSpace(:) == 1);

end