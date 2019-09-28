function ppicaps_toyExample
%
% PPICAPS_TOYEXAMPLE - This function is meant to illustrate how the 
% PPI-CAPs time-resolved analysis works. For more information on PPI-CAPs
% please check (Freitas et al., 2019).

addpath(genpath([pwd '/dependencies']));


% -------------------------------------------------------------------------
% SET UP SIMULATED EXPERIMENTAL DATA
% -------------------------------------------------------------------------

% Let us suppose that we have some experimental task-based data. We would
% like to find out how the relationship between a particular seed activity
% and the organisation of the rest of the brain evolve over time. Are there
% patterns whose correlation with the seed is higher during one specific
% task? And are there others whose relationship with the seed activity is
% constant independently of the task being performed? At what specific
% times do these brain patterns actually happen, throughout the experiment?

% To answer those questions as per the PPI-CAPs approach, we start by 
% selecting time points where the seed  is most highly active or deactive. 
% This selection can be done based on the percentage of frames you would
% like to keep in the analysis, or based on a threshold of your choice ? as
% we show on the paper, the method is pretty robust to a wide range of
% threshold choises. Since the seed activity must be orthogonal to the task
% timecourse, this will mean that we end up selecting a similar number of 
% frames where the seed is active or deactive for each of the tasks.

% Below, are our true labels for the sign of the seed activity, the task 
% timecourse and the ppi (as an element-wise multiplication of the seed and
% task signs), on the suprathreshold frames. That is, we have 24 frames that
% were selected at moments where the seed was highly active or deactive.
seedlabels = [1 -1 1 -1 1 -1  1 -1  1 -1  1 -1 1 -1 1 -1 1 -1  1 -1  1 -1  1 -1];
tasklabels = [1  1 1  1 1  1 -1 -1 -1 -1 -1 -1 1  1 1  1 1  1 -1 -1 -1 -1 -1 -1];
ppilabels  = [1 -1 1 -1 1 -1 -1  1 -1  1 -1  1 1 -1 1 -1 1 -1 -1  1 -1  1 -1  1];


% -------------------------------------------------------------------------
% SET UP SIMULATED GROUND TRUTH
% -------------------------------------------------------------------------

% Let us assume that there are three brain patterns, composed of 5 voxels 
% each, as such: 
ppicap1 = [-0.7 -0.75 -0.85 -0.9 -0.95];
ppicap2 = [0.8 0.85 0.9 -0.8  -0.85];
ppicap3 = [-0.7 0.9 -0.75 0.95 -0.7];

% Let's plot them and check them out.
figure; 
colors = cbrewer('div', 'Spectral', 50, 'pchip'); caxis([-1.5 1.5]);
imagesc([ppicap1(:),ppicap2(:), ppicap3(:)]);colormap(colors);caxis([-2 2]);
set(gca, 'xtick', 1:3, 'xticklabels',{'Seed', 'Task', 'PPI'},'XTickLabelRotation',90,'FontSize',10);
set(gcf, 'color', 'w', 'pos', [-1562 292 190 473]);colorbar; 


% At different time points, each of these patterns may be expressed with 
% one of two polarities: positive or negative. That is, when ppicap1 is 
% expressed with a negative polarity, every voxel has the sign of its 
% intensity flipped as compared to its positive polarity. In practical 
% terms, this would mean that regions that are active in one polarity are 
% de-active in the other polarity, and vice-versa. 

% Moreover, the polarity of each pattern may be random, or it may correlate 
% with the  activity of a seed region, the timecourse of the task, or an 
% interaction of the two, that is: it correlates more with the seed during 
% one of the tasks (or it correlates more with the task when the seed is 
% either active or deactive). Let's call this interaction a PPI 
% (psychophysiological interaction). And we will represent this as shown 
% below:

sp = ppicap1(:); %seedpattern (we will force this to correlate with the seed)
tp = ppicap2(:); %taskpattern (we will force this to correlate with the task)
pp = ppicap3(:); %ppipattern (we will force this to correlate with the ppi)
patternsTimeCourse = [tp,  tp, sp, -sp, pp, -pp,...
                     -tp, -tp, sp, -sp, -pp, pp,...
                      tp,  tp, sp, -sp, pp, -pp,...
                     -tp, -tp, sp, -sp, -pp, pp];

% So let us visualise these patterns. Each one has been labelled with the 
% effect (seed / task / ppi, in positive or negative polarity) we have 
% forced them to correlate with, for reference purposes.
figure; imagesc(patternsTimeCourse);colormap(colors);
patternsTimeCourseLabels = {'task', 'task', 'seed', '-seed', 'ppi', '-ppi',...
                          '-task', '-task', 'seed', '-seed', '-ppi', 'ppi',...
                            'task', 'task', 'seed', '-seed', 'ppi', '-ppi',...
                          '-task', '-task', 'seed', '-seed', '-ppi', 'ppi'};
set(gca, 'xtick', 1:24, 'xticklabels',patternsTimeCourseLabels,'XTickLabelRotation',90,'FontSize',10);set(gcf, 'color', 'w');caxis([-2 2]);



% -------------------------------------------------------------------------
% STATIONARY ANALYSIS (SiMAP)
% -------------------------------------------------------------------------

% This is what a stationary interaction map would look like, when we
% multiply frames by the sign of the corresponding PPI values, and then add
% those together:
figure; imagesc(sum(patternsTimeCourse.*repmat(ppilabels, 5, 1),2));colormap(colors);set(gcf, 'color', 'w');caxis([-16 16])

% Note that, in this case, since we only had one pattern that had a PPI
% effect, the stationary analysis was able to recover it perfectly.
% However, if we had two or more different patterns that correlated with 
% the PPI time course, the stationary map would return an average of those,
% therefore depicting a pattern that never actually happens in the brain,
% at any point in time!!! This is why we need to proceed with a dynamic
% analysis.

% -------------------------------------------------------------------------
% DYNAMIC ANALYSIS ATTEMPTS
% -------------------------------------------------------------------------


% Now let us try to recover the underlying patterns by clustering the 
% suprathreshold frames using regular k-means, where k = 3 (since we know
% that this is the exact number of underlying patterns):
[IDX, C] = kmeans(patternsTimeCourse',3,'dist', 'cosine', 'replicates', 100 );
figure;imagesc(C');colormap(colors); title('K-means (K = 3)');set(gcf, 'color', 'w');caxis([-1 1]);

% Note that k-means with k = 3 returns patterns that don't even originally
% exist in our dataset, so this would not be a viable option.


% Now, let us try again to recover the underlying patterns by clustering 
% the frames using k-means with k = 6:
[IDX, C] = kmeans(patternsTimeCourse',6,'dist', 'cosine', 'replicates', 100 );
figure;imagesc(C');colormap(colors);title('K-means (K = 6)');set(gcf, 'color', 'w');

% Note that now we have indeed recover all the different occurences of our 
% three patterns. However, when a pattern reappears with a different polarity,
% the two occurences are assigned to different centroids! This would not 
% make much sense as the pattern is the same up to a sign flip. More
% importantly, we would need to create and additional metric to
% automatically match centroids that are opposite polarities of the same
% pattern. 


% -------------------------------------------------------------------------
% DYNAMIC ANALYSIS WITH PPI-CAPS
% -------------------------------------------------------------------------

% So let us instead use PPI-CAPs and cluster using the k-means based on
% modulo-pi cosine distance, to retrieve the flipping of each frame at
% clustering time.
k = 3;
[PPI_CAPs, idx] = kmeanspp(patternsTimeCourse, k, 'cosine', 100); 

for i = 1:k
    ppi_cap_vols{i} = patternsTimeCourse(:,find(PPI_CAPs == i)); %instead of frames use test = patternsTimeCourse.*repmat(sign(ppi(idxSupraThresh)), 9,1);
    ppi_cap_pols{i} = idx(PPI_CAPs == i); ppi_cap_pols{i}(ppi_cap_pols{i}==2) =-1; 
    ppi_cap(:,i) = sum(ppi_cap_vols{i}.*repmat(ppi_cap_pols{i}, 5,1),2);
end

% Now note how we manage to recover the exact patterns we started out with,
% up to a sign flip. Inportantly, the sign itself does not imply that the
% pattern was negative, but that at each given timepoint, each voxel must
% be multiplied by the corresponding sign of the PPI-CAP so we can trully
% know if it is active or deactive at that point in time.
figure;imagesc(ppi_cap); colormap(colors); 
title('K-means with modulo-pi dist.');caxis([-16 16]);



% -------------------------------------------------------------------------
% STATISTICAL ANALYSIS
% -------------------------------------------------------------------------

% So now let us suppose we did not know what the group truth was, and let 
% us use confusion matrices and permutation testing to check whether each 
% PPI-CAP's frames correlate with any of the effects we are trying to 
% analyse.

 
seedSigns = seedlabels;
taskSigns = tasklabels;

maxPerm = 1000; % Maximum number of random permutations
figure;

% Calculate confusion matrices for each PPI-CAP.
for i = 1:3
    thisCapSeed = seedSigns(PPI_CAPs==i);
    thisCapTask = taskSigns(PPI_CAPs==i);
    thisCapPPI = thisCapSeed .* thisCapTask;
    thisCapClusteringFlip = idx(PPI_CAPs==i);
    
    % Seed
    subplot(3,3,3*i-2);
    [seedConfusionMatrix{i}, seed_fkix{i}] =  calc_confusionMatrix(thisCapSeed, thisCapClusteringFlip);
    imagesc(seedConfusionMatrix{i}); colormap(colors);
    
    % Task
    subplot(3,3,3*i-1);
    [taskConfusionMatrix{i},  task_fkix{i}] = calc_confusionMatrix(thisCapTask, thisCapClusteringFlip);
    imagesc(taskConfusionMatrix{i});colormap(colors);
     
    % PPI
    subplot(3,3,3*i);
    [ppiConfusionMatrix{i}, ppi_fkix{i}]  = calc_confusionMatrix(thisCapPPI, thisCapClusteringFlip);
    imagesc(ppiConfusionMatrix{i});colormap(colors);
    
    
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


% Now let us plot the distributions of fk-indeces obtained through random 
% permutations of the (seed; task; and ppi) labels from the frames that
% make up each PPI-CAP. Then, we plot a line for the real data-driven 
% fk-index and see where it stands!
figure;
for i = 1:3
    
    % Seed
    ax2 = subplot(3,3,3*i-2);
    h2 = histogram(ax2,seed_fkix_rnd(i,:));title(ax2, 'Null Seed Distribution');
    h2.BinMethod = 'Sturges'; h2.Normalization = 'countdensity'; h2.FaceColor = colors(i,:); h2.NumBins = 5;
    line(ax2, [seed_fkix{i}, seed_fkix{i}], get(ax2, 'ylim'), 'LineWidth', 2, 'Color', 'r');
    pSeed = min(((length(find(seed_fkix_rnd(i,:)<seed_fkix{i}))+1)/(maxPerm+1)),((length(find(seed_fkix_rnd(i,:)>seed_fkix{i}))+1)/(maxPerm+1)));
    text(ax2, 0.7,0.9,['p = ' num2str(pSeed)],'Units','normalized');
    
    % Task
    ax4 = subplot(3,3,3*i-1);
    h2 = histogram(ax4,task_fkix_rnd(i,:));title(ax4, 'Null Task Distribution');
    h2.BinMethod = 'Sturges'; h2.Normalization = 'countdensity'; h2.FaceColor = colors(i,:);h2.NumBins = 5;
    line(ax4, [task_fkix{i}, task_fkix{i}], get(ax4, 'ylim'), 'LineWidth', 2, 'Color', 'r');
    pTask = min(((length(find(task_fkix_rnd(i,:)< task_fkix{i}))+1)/(maxPerm+1)),((length(find(task_fkix_rnd(i,:)>task_fkix{i}))+1)/(maxPerm+1)));
    text(ax4, 0.7,0.9,['p = ' num2str(pTask)],'Units','normalized');
    
    
    % PPI
    ax3 = subplot(3,3,3*i);
    h2 = histogram(ax3, ppi_fkix_rnd(i,:));title(ax3, 'Null PPI Distribution');
    h2.BinMethod = 'Sturges'; h2.Normalization = 'countdensity'; h2.FaceColor = colors(i,:);h2.NumBins = 5;
    line(ax3, [ppi_fkix{i}, ppi_fkix{i}], get(ax3, 'ylim'), 'LineWidth', 2, 'Color', 'r');
    pPPI = min(((length(find(ppi_fkix_rnd(i,:)< ppi_fkix{i}))+1)/(maxPerm+1)),((length(find(ppi_fkix_rnd(i,:)>ppi_fkix{i}))+1)/(maxPerm+1)));
    text(ax3, 0.7,0.9,['p = ' num2str(pPPI)],'Units','normalized');
    
end



end


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

confusionMatrix = [sum((thisCapsEffectLabels==1).*(thisCapClusteringFlip ==1)), ...
    sum((thisCapsEffectLabels==1).*(thisCapClusteringFlip ==2)); ...
    sum((thisCapsEffectLabels==-1).*(thisCapClusteringFlip ==1)), ...
    sum((thisCapsEffectLabels==-1).*(thisCapClusteringFlip ==2))];

fk_index = det(confusionMatrix)/((numel(thisCapsEffectLabels)/2)^2);

end

