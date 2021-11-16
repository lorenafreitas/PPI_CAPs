%% This function performs consensus clustering over a range of K values
% The goal is to provide a measure of how good each value of K is
% 
% Inputs:
% - X is the data matrix (n_DP x n_DIM)
% - K_range is the range of K values to examine
% - Subsample_type defines how subsampling is done: across items (data
% points) if 'items', and across dimensions if 'dims', across subjects (all
% items of one subject) if 'subjects'
% - Subsample_fraction is the fraction of the original data points, or
% dimensions, to keep for a given fold
% - n_folds is the number of folds over which to run
% - DistType: 
% - subject_labels is the assignment of every item in X to a subject, only
% needed for Subsample_type='subjects'
%
% Outputs:
%   consensusOutput - struct containing consensus infromation
%       .K_range
%       .Consensus_ordered - Consensus matrix for every fold (n_items x
%       n_items x nK)
%       .clusters - struct containing cluster centers, IDX etc for every fold
%
%
% 15.11.2017 - Lorena:
   %    created based on Daniela's version, added use of kmeansCustom
%
%
% 14.6.2017 - Daniela:
%       created based on Thomas' version, added subject subsampling,
%       changed outputs to structs
%
% 10.9.2021 - Updated & modified by: Serafeim Loukas on Sep 10, 2021
%       customization: see lines 205 & 249


function [consensusResults] = ConsensusClustering(X,subject_labels,param,outDir_cons, case_name)
    
    K_range=param.K;
    n_folds=param.cons_n_folds;
    DistType=param.DistType;
    Subsample_fraction=param.Subsample_fraction;
    Subsample_type=param.Subsample_type;
   
    if isfield(param,'MaxIter')
        MaxIter=param.MaxIter;
    else
        MaxIter=100;
    end
    
    if isfield(param,'KmeansMethod')
        KmeansMethod=param.KmeansMethod;
    else
        KmeansMethod='kmeanspp';
    end
    
    
fprintf(['\nRunning consensus clustering using ' KmeansMethod ' method \n' ...
    '--------------------------------------------------------\n\n']);
    
    clear param
    
    % Number of data points
    n_items = size(X,1);
    
    % Number of dimensions
    n_dims = size(X,2);
    
%     Consensus = sparse(n_items,n_items);
%     Consensus_ordered = sparse(n_items,n_items,length(K_range));
    
    
    
    % Loop over all K values to assess
    for k = 1:length(K_range)
        
        disp(['Running consensus clustering for K = ',num2str(K_range(k)),'...']);
        
        % Connectivity matrix that will contain 0s or 1s depending on whether
        % elements are clustered together or not
%         M = sparse(n_items,n_items,n_folds);
%         I = sparse(n_items,n_items,n_folds);
        M_sum=zeros(n_items,n_items);
        I_sum=zeros(n_items,n_items);
        
%         disp('before h loop');
        
        % Loops over all the folds to perform clustering for
        for h = 1:n_folds
            disp(['Fold ' num2str(h) ':'])
            switch Subsample_type
                case 'items'
                    
                    % Number of items to subsample
                    n_items_ss = floor(Subsample_fraction*n_items);
                    
                    % Does the subsampling
                    [X_ss,tmp_ss] = datasample(X,n_items_ss,1,'Replace',false);
                    
                    % Vector
                    I_vec = zeros(n_items,1);
                    I_vec(tmp_ss) = 1;  
                    
                    I_vec_s=sparse(I_vec);
                    
                    %Constructs the indicator matrix
                    I_sum=I_sum+I_vec_s*I_vec_s';
                    
                    
                case 'dims'
                    
                    % Number of dimensions to subsample
                    n_dims_ss = floor(Subsample_fraction*n_dims);
                    
                    % Does the subsampling
                    [X_ss,tmp_ss] = datasample(X,n_dims_ss,2,'Replace',false);
                    
                    % Constructs the indicator matrix
                    I(:,:,h) = ones(n_items,n_items);
                    
                case 'subjects'
                    
                    if ~exist('subject_labels','var')
                        error('subject labels needed for subjects subsampling!')
                    end
                    
                    subject_list=unique(subject_labels);
                    n_subjects=length(subject_list);
                    
                    % Number of items to subsample
                    n_subjects_ss = floor(Subsample_fraction*n_subjects);
                    
                    % Does the subject subsampling
                    disp('Subject subsampling...'); tic
                    [subjects_ss,~] = datasample(subject_list,n_subjects_ss,1,'Replace',false);
                    
                    % Vector
                    I_vec = zeros(n_items,1);
                    for iS=1:n_subjects_ss
                        I_vec(subject_labels==subjects_ss(iS)) = 1;
                    end
                    
                    I_vec_s=sparse(I_vec);
                    
                    %Constructs the indicator matrix
                    I_sum=I_sum+I_vec_s*I_vec_s';
                    toc
                    
                    % subsampled data
                    disp('Data subsampling...');tic
                    X_ss=X(I_vec>0,:);
                    toc
                otherwise
                    errordlg('PROBLEM IN TYPE OF SUBSAMPLING');
            end
            
            % Does the clustering (for now, only with k-means), so that IDX
            % contains the indices for each datapoint
            disp('Clustering ...');tic
            %IDX = kmeans(X_ss,K_range(k),'Distance',DistType,'Replicates',10,'Display','final','MaxIter',MaxIter);
            %IDX   = kmeansCustom(X_ss', K_range(k), 10);
            switch KmeansMethod
                case 'kmeanspp'
                         IDX   = kmeanspp(X_ss', K_range(k), DistType, 10); IDX = IDX';
                case 'kmeansmatlab'
                         IDX   = kmeans(X_ss, K_range(k),'Distance', DistType, 'Replicates', 10); IDX = IDX';
            end
            toc
            
            % just to check whether it would be faster to save every
            % iteration right away, instead of using the working memory
            if ~exist([outDir_cons '/tmp/'],'dir'); mkdir([outDir_cons '/tmp']);end
            save([outDir_cons '/tmp/IDX_' num2str(h)],'IDX');
            save([outDir_cons '/tmp/I_vec_' num2str(h)],'I_vec');
            %save([outDir_cons '/tmp/subjects_ss_' num2str(h)],'subjects_ss');
            
            clear X_ss
            
            % Builds the connectivity matrix
            disp('Buiding connectivity matrix M ...');tic
            M_sum=M_sum+Build_Connectivity_Matrix(IDX,find(I_vec>0),Subsample_type,n_items);
            toc
            
            
        end
        
        
        disp('Computing Consensus Matrix ...');tic
        % compute sum for sparse cell arrays
%         M_sum=sparse(n_items,n_items);
%         I_sum=sparse(n_items,n_items);
%         for h = 1:n_folds
%             M_sum=M_sum+M{h};
%             I_sum=I_sum+I{h};
%         end
        
        % Constructs the consensus matrix for the considered K
        Consensus = M_sum./I_sum;
        if any(I_sum(:)==0)
            warning([num2str(nnz(isnan(Consensus(:)))) ' (' ...
                num2str(nnz(isnan(Consensus(:)))/length(Consensus(:))) ...
                '%) items have not been selected during subsampling, you should increase the number of folds!']);
            Consensus(isnan(Consensus))=0;
            % MODIFIED by SERAFEIM: diag should have 1s and not a mix of 1s and 0s due to NaN masking
            Consensus(logical(eye(size(Consensus,1)))) = 1;
        end
        toc
        
        disp('Ordering consensus matrix ...'); tic
        tree = linkage(squeeze(1-Consensus),'average');

        % Leaf ordering to create a nicely looking matrix
        leafOrder = optimalleaforder(tree,squeeze(1-Consensus));
        
        % Ordered consensus matrix
        Consensus_ordered = Consensus(leafOrder,leafOrder);
        toc
        
        % load IDX and I_vec data
        disp('Loading and concatenating data...');tic
        for h = 1:n_folds
            load([outDir_cons '/tmp/IDX_' num2str(h) '.mat']);
            IDX_all{h}=IDX;
            load([outDir_cons '/tmp/I_vec_' num2str(h) '.mat']);
            I_vec_all(:,h)=I_vec;
           % load([outDir_cons '/tmp/subjects_ss_' num2str(h) '.mat']);
            %subjects_ss_all(:,h)=subjects_ss;
        end
        IDX=IDX_all;
        I_vec=I_vec_all;
        %subjects_ss=subjects_ss_all;
        toc
        
        disp('Saving data in result struct ...');tic
        consensusResults{k}.K=K_range(k);
        %consensusResults{k}.subjects_ss=subjects_ss;
        consensusResults{k}.I_vec=I_vec;
        consensusResults{k}.IDX=IDX;
%         consensusOutput{k}.I=I;
%         consensusOutput{k}.M=M;
%         consensusOutput{k}.I_sum=I_sum;
%         consensusOutput{k}.M_sum=M_sum;
%         consensusOutput{k}.Consensus=Consensus;
        consensusResults{k}.Consensus_ordered=Consensus_ordered;
        consensusResults{k}.leafOrder=leafOrder;
        toc
        
        %save(fullfile(outDir_cons,'consensusResults'),'consensusResults','-v7.3')
        save(fullfile(outDir_cons,['consensusResults' '_thres_' case_name 'perc']),'consensusResults','-v7.3')
        
        clearvars leafOrder M I M_sum I_sum Consensus Consensus_ordered
    end
end