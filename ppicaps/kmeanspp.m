function [L,Ix,C,dis,sD] = kmeanspp(X,k, dist,varargin)
%   KMEANSPP ? Cluster multivariate data using the k-means++ algorithm.
%   [L,Ix,C,dst,sD] = kmeans(X,k) takes an n-by-p matrix X as input containing p observations of
%   n-dimensional data and produces a 1-by-size(X,2) vector L with one class
%   label per column in X, a size(X,1)-by-k matrix C containing the
%   centers corresponding to each class,  a vector Ix indicating the
%   sign of the time point (1 = positive; 2 = negative) a vector dis 
%   containing the average distance from each centroid to their corresponding 
%   points and a vector sD indictating the standard deviation of each
%   cluster
%
%   Author: Lorena Freitas (Lorena.Freitas@epfl.ch)
%            
%   %   Last checked: September 2019
%
%   'varargin'
%            #replicates
%            or C - matrix of centroids
%
%   References:
%   [1] J. B. MacQueen, "Some Methods for Classification and Analysis of
%       MultiVariate Observations", in Proc. of the fifth Berkeley
%       Symposium on Mathematical Statistics and Probability, L. M. L. Cam
%       and J. Neyman, eds., vol. 1, UC Press, 1967, pp. 281-297.
%   [2] D. Arthur and S. Vassilvitskii, "k-means++: The Advantages of
%       Careful Seeding", Technical Report 2006-13, Stanford InfoLab, 2006.


% reseed random number generator
%rng shuffle

% Default value of replicates: only used when centroids are passed as an
% argument, which means we only need to assign each frame to a centroid,
% and therefore there's no need to run it more than once.
replicates = 1;

switch dist
    case 'sqEuclidean'
        scoremethod    = @sqEuclideanScore;
        distancemethod = @sqEuclideanDistance;
    case 'cosine'
        scoremethod    = @CosineScore;
        distancemethod = @CosineDistance;
        
        % Normalise vectors to unit norm
        X = X./repmat(sqrt(sum(X.*X)),size(X,1),1); 
end

if nargin > 3
    
    % If 3rd argument is a matrix of centroids, assign frames to centroids
    if length(varargin{1}) > 1
        C         = varargin{1};  
        [M1,Lpos] = min(scoremethod(C,X),[],1);
        [M2,Lneg] = min(scoremethod(-C,X),[],1);
        [~,Ix]    = min([M1;M2]);
        L         = ones(1,length(Lpos)); L(Ix==1) = Lpos(Ix==1); L(Ix==2) = Lneg(Ix==2);
        M         = ones(1,length(M1));  M(Ix==1) = M1(Ix==1); M(Ix==2) = M2(Ix==2);
   
    % if 3rd argument is the number of Ks
    else
        replicates          = varargin{1};
        bestSolution.optVar = Inf; %realmin;
        
        for repl=1:replicates
            
            L = [];
            L1 = 0;
            
            while length(unique(L)) ~= k
                
                % The k-means++ initialization. The code below is inspired
                % by Laurent Sorber's code at MATLAB Central File Exchange
                % from  2013-02-08.
                C = X(:,1+round(rand*(size(X,2)-1)));
                L = ones(1,size(X,2));
                for i = 2:k
                    D = X-C(:,L); D2 = X+C(:,L); % D2 is distance to negative centroid;
                    D = min([cumsum(sqrt(dot(D,D,1)));cumsum(sqrt(dot(D2,D2,1)))]);
                    if D(end) == 0, C(:,i:k) = X(:,ones(1,k-i+1)); return; end
                    C(:,i) = X(:,find(rand < D/D(end),1));
                    
                    [M1,Lpos] = min(scoremethod(X,C));
                    [M2,Lneg] = min(scoremethod(X,-C));
                    [~,Ix]    = min([M1;M2]);
                    L         = ones(1,length(Lpos)); L(Ix==1) = Lpos(Ix==1); L(Ix==2) = Lneg(Ix==2);
                end
                
                % The actual k-means algorithm.
                while any(L ~= L1)
                    L1 = L;
                    
                    % Update clusters
                    for i = 1:k, l = L==i; 
                        scf         = l.*Ix; % signed cluster's frames (1 = pos; 2 = neg)
                        C(:,i)      = (sum(X(:,scf==1),2) - sum(X(:,scf==2),2))/sum(l);
                         
                    end
                    
                    [M1,Lpos] = min(scoremethod(X,C),[],1);
                    [M2,Lneg] = min(scoremethod(X,-C),[],1);
                    [~,Ix]    = min([M1;M2]);
                    L         = ones(1,length(Lpos)); L(Ix==1) = Lpos(Ix==1); L(Ix==2) = Lneg(Ix==2);
                end
                
                M          = ones(1,length(M1));  M(Ix==1) = M1(Ix==1); M(Ix==2) = M2(Ix==2);
                
                % Calculate sumD
                optVar     = sum(M);
                
                % Is this the best solution so far? If so, remember this one.
                if optVar < bestSolution.optVar
                    bestSolution.optVar = optVar;
                    bestSolution.L    = L;
                    bestSolution.C    = C;
                    bestSolution.Ix   = Ix;
                end
            end
        end
        % Return best solution
        L      = bestSolution.L;
        C      = bestSolution.C;
        Ix     = bestSolution.Ix;
         
    end
    
    % calculate actual distances
    dstpos = distancemethod(X,C); 
    dstneg = distancemethod(X,-C);
    dst    = ones(size(dstpos)); dst(:,Ix==1) = dstpos(:,Ix==1); dst(:,Ix==2) = dstneg(:,Ix==2);

    dis = nan(1,k);
    sD = nan(1,k);
    for i = 1:k
        dis(i) = sum(dst(i,L==i));
        sD(i)  = std(dst(i,L==i));
    end
end
end


%% ------------------------------------------------------------------------
%  DISTANCE FUNCTIONS
%  ------------------------------------------------------------------------
%  These functions have been optimised to identify the nearest cluster in a
%  vectorised way.
%  ------------------------------------------------------------------------
function dist = sqEuclideanDistance(X,C)

dist = bsxfun(@plus,dot(X,X,1), sqEuclideanScore(X,C));

end

function dist = CosineDistance(X,C)

dist = 1 - bsxfun(@rdivide,(X'*C),(sqrt(dot(X,X,1))'*sqrt(dot(C,C,1))))';

end

%% ------------------------------------------------------------------------
%  DISTANCE SCORE FUNCTIONS
%  ------------------------------------------------------------------------
%  These functions have been optimised to identify the nearest cluster in a
%  vectorised way using the minimum computation possible. Their output
%  should be thought of as a score rather than a distance  metric.
%  ------------------------------------------------------------------------
function score = sqEuclideanScore(X,C)

score = bsxfun(@minus, dot(C,C,1), 2*X'*C)';

end

function score = CosineScore(X,C)
% multiply the score by -1 so we optimise for the minimum
% value  as opposed to maximum (i.e.: distance, not similarity)
 
score = bsxfun(@rdivide, (C'*X)', sqrt(dot(C,C,1)))'*-1;

end
