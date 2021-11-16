%% function to compute the quality indices for the comparison of different k
%
% 15.6.2017 - Daniela Z?ller: created based on function by Thomas
%               adapted for large sparse matrices structures

function [CDF,AUC,Delta] = ComputeClusteringQuality(consensusResults,K_range, plotFlag)


if nargin <3
    plotFlag = 0;
end
% Number of K values to check
nK = length(K_range);

if iscell(consensusResults)
    n_items = size(consensusResults{1}.Consensus_ordered,1);
else
    n_items = size(consensusResults,1);
end

% Creates the CDF range
c = 0:0.01:1;

colors = cbrewer('qual', 'Set3', nK);
set(groot,'defaultAxesColorOrder',colors);


figure;
hold on;


for k = 1:nK
    disp(num2str(K_range(k)))
    if iscell(consensusResults)
        Consensus_k=consensusResults{k}.Consensus_ordered;
    else
        Consensus_k=squeeze(consensusResults(:,:,k));
    end
    
    
    % Sorted consensus entries
    Cons_val = jUpperTriMatToVec(Consensus_k);
    Cons_val = sort(Cons_val,'ascend');
    
    % Computation of CDF
    for i = 1:length(c)
        %             if ~mod(i,10);disp(num2str(i));end
        CDF(k,i) = nnz(Cons_val <= c(i));
    end
    CDF(k,:)=CDF(k,:)./length(Cons_val);
    
    
    % vectorized computation of AUC
    %         Cons_val_diff=diff(Cons_val);
    %         CDF_interp=interp1q(c',CDF(k,:)',Cons_val);
    AUC(k)=diff(c)*CDF(k,2:end)';
    
    
    pl(k)=plot(c,CDF(k,:),'linewidth',2);
    legend(pl,[num2str(K_range')],'location','southeast'); ylabel('CDF'); xlabel('Consensus Index Value');set(gcf,'color','white');grid on;
    % Computation of the AUC
    %         AUC(k) = 0;
    %         for i = 2:length(Cons_val)%(n_items*(n_items-1)/2)
    %             AUC(k) = AUC(k) + (Cons_val(i)-Cons_val(i-1))* interp1q(c',CDF(k,:)',Cons_val(i));
    %         end
    
    %         clear Cons_val CDF_interp Cons_val_diff
end





% Computation of Delta AUC
max_AUC = AUC(1);
Delta(1) = AUC(1);
for k = 2:nK
    Delta_percent(k) = (AUC(k) - max_AUC)/max_AUC;
    Delta(k) = AUC(k) - max_AUC;
    if AUC(k)>max_AUC
        max_AUC=AUC(k);
    end
end



%% plotting results
%     figure;
%     hold on;
%     for k = 1:nK
%         plot(c,CDF(k,:));
%     end
%     legend(cellstr(num2str(K_range')));
%title('CDF');
if plotFlag
    figure;
    plot(K_range(1:length(AUC)),AUC,'linewidth',2);
    title('AUC');
    
    figure;
    plot(K_range(2:length(Delta)),Delta(2:end),'linewidth',2);
    title('Delta AUC');
    
    figure;
    plot(K_range(2:length(Delta_percent)),Delta_percent(2:end),'linewidth',2);
    title('Delta AUC (percentage)');
    
end
end