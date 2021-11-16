% This function calculates the confusion matrix for each PPI-CAP. The idea
% is that, if a coactivation pattern has a specific effect (seed / task / 
% ppi), the polarity of its composing frames as returned by the clustering
% method should correlate with the sign of that effect. So, for each effect,
% we count the number of frames within that PPI-CAP that are positive when
% the timecourse of that effect is positive, and vice-versa, obtaining a
% confusion matrix. If there is a high correlation between the two, the
% determinant of the matrix will have a high value. The fk-index is thus
% the determinant of the matrix normalised by the total number of frames 
% within the PPI-CAPs.
function [confusionMatrix, fk_index] =  calc_confusionMatrix(thisCapsEffectLabels, thisCapClusteringFlip)

confusionMatrix = [sum((thisCapsEffectLabels==1) .*(thisCapClusteringFlip ==1)), ...
                   sum((thisCapsEffectLabels==1) .*(thisCapClusteringFlip ==2)); ...
                   sum((thisCapsEffectLabels==-1).*(thisCapClusteringFlip ==1)), ...
                   sum((thisCapsEffectLabels==-1).*(thisCapClusteringFlip ==2))];

fk_index = det(confusionMatrix)/((numel(thisCapsEffectLabels)/2)^2);

end