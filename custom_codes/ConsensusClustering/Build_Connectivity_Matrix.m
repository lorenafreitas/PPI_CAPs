%% This function constructs a connectivity matrix from the indices of
% cluster assignment
function [M] = Build_Connectivity_Matrix(IDX,tmp_ss,type,n_items)

    IDX_full = zeros(n_items,1);

    switch type
        case 'items'
            IDX_full(tmp_ss) = IDX;
        case 'subjects'
            IDX_full(tmp_ss) = IDX;
        case 'dims'
            IDX_full = IDX;
    end
    
    M = sparse(n_items,n_items);
    for iK=1:length(unique(IDX)) % for every cluster
        IDX_tmp=sparse(n_items,1);
        IDX_tmp(IDX_full==iK)=1;
        M=M+IDX_tmp*IDX_tmp'; % matrix with all elements that are in cluster i
    end
    
    
%     IDX_full=sparse(IDX_full);
%     M=IDX_full*IDX_full';
    
%     M = zeros(n_items,n_items);
% 
%     for i = 1:length(IDX_full)
%         for j = 1:length(IDX_full)
%             if (IDX_full(i) == IDX_full(j)) && (IDX_full(i) > 0)
%                 M(i,j) = 1;
%             end
%         end
%     end
end