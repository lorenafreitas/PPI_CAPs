function consensusClusteringPlotResults(consensusResults, resDir, case_name)
% -------------------------------------------
%    Updated & modified by Serafeim Loukas on July 19, 2021
% -------------------------------------------
warning('off')

%% Clustering Quality Assessment
% -------------------------------
K_range           = consensusResults{1}.K:consensusResults{end}.K;
[CDF2,AUC2,Delta2] = ComputeClusteringQuality(consensusResults,K_range);


%% Cumulative Distribution Function (CDF)
% -----------------------------------------------------
% (Monti et al., Machine Learning, 2003)
% -----------------------------------------------------
set(0,'defaultAxesFontSize',14)
nK=length(K_range);
c = 0:0.01:1;
colors = cbrewer('seq', 'GnBu', nK+5); colors=colors(6:end,:);
colors2=cbrewer('qual', 'Set1',nK); colors(9,:)=colors2(1,:);
set(groot,'defaultAxesColorOrder',colors);
figure('position',[440   378   515   420]);
hold on;

% Plot one line per value of k
for k = 1:nK
    pl(k)=plot(c,CDF2(k,:),'linewidth',2,'color',colors(k,:));
    legend(pl,[num2str(K_range')],'location','southeast');
end

% Edit plot's appearance
title('CDF');ylabel('CDF'); xlabel('Consensus Index Value');set(gcf,'color','white');grid on;


%% Proportion of Ambiguous Clustering (PAC)
% -----------------------------------------------------
% (Senbabaoglu et al.,PLOS Computational Biology, 2014)
% -----------------------------------------------------
figure('position',[440 378 216 314]);
for k = 1:nK
    pac(k+1) = CDF2(k,90) - CDF2(k,10); 
end
%plot(2:k+1, pac(2:end),'linewidth',2,'color',colors(1,:)); grid on; hold on;
plot(K_range, pac(2:end),'linewidth',2,'color',colors(1,:)); grid on; % FIXED: by Serafeim Loukas, July 16, 2021
ylabel('PAC (Proportion of Ambiguous Clustering)'); xlabel('Number of Clusters K')


%% Area Under the Curve (AUC) of the CDF
% -----------------------------------------------------
% (Monti et al., Machine Learning, 2003)
% -----------------------------------------------------
figure('position',[440   378   216   314]);
hold on;
plot(K_range,AUC2,'linewidth',2,'color',colors(1,:)); grid on;
set(gcf,'color','white');
ylabel('AUC'); xlabel('Number of Clusters K');

%% Delta AUC (Change in AUC from one value to the next)
% -----------------------------------------------------
% (Monti et al., Machine Learning, 2003)
% -----------------------------------------------------
figure('position',[440 378   216    314]);
hold on;
a1 = gca;
plot(K_range(1:end),Delta2(1:end),'linewidth',2,'color',colors(1,:));
grid on;
ylabel('\Delta AUC'); xlabel('Number of Clusters K')
set(gcf,'color','white');
%print('Delta','-depsc2','-painters');

%% Consensus Matrix, for Visual Inspection
% -----------------------------------------------------
% (Monti et al., Machine Learning, 2003)
% -----------------------------------------------------
% Heatmap of consensus matrix after sorting using hierarchical clustering
disp('Heatmap of consensus matrix after sorting using hierarchical clustering')
figure('position', [  237   367   883   402]);
thisCmap = cbrewer('div', 'RdYlBu', 50);
for iK=1:length(consensusResults);
    subplot(2,4,iK);
    im = imagesc(consensusResults{iK}.Consensus_ordered);
    colormap(flipud(thisCmap));
    im.AlphaData = 0.9;
    title(['k=' num2str(consensusResults{iK}.K)]);
end
box off; set(gcf,'color','white');


%% Save results as .png figures
resDir=fullfile(resDir,'plots');
mkdir(resDir);
h = get(0,'children');
for i=1:length(h)
    print(h(i), fullfile(resDir,['figure_' num2str(i) '_thres_' case_name 'perc']), '-dpng', '-r300', '-painters');
end

end