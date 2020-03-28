% Function based on the relevant section of SPM8's spm_peb_ppi.m 
% function for devonvolution, usually used in the context of PPI analyses.
%
% Last checked: September 2019

function deconvData = deconvolveSubj(subj)


% Load SPM
%--------------------------------------------------------------------------
analyze_dir = [subj.dataDir 'GLM_SciFun/'];
SPM = spm_select('FPListRec',analyze_dir,'SPM.mat$'); load(SPM);
tempVOI = load([subj.dataDir '/VOI_PCC_1.mat']);
xY = tempVOI.xY;
VI =  SPM.xY.VY;

% Setup variables
%--------------------------------------------------------------------------
RT      = SPM.xY.RT;
dt      = SPM.xBF.dt;
NT      = round(RT/dt);
N       = length(xY(1).u);
k       = 1:NT:N*NT;
Sess    = SPM.Sess(1);

% Create basis functions and hrf in scan time and microtime
%--------------------------------------------------------------------------
hrf = spm_hrf(dt); % 1 corresponds to the microtime resolution I want


% Create convolved explanatory {Hxb} variables in scan time
%--------------------------------------------------------------------------
xb  = spm_dctmtx(N*NT + 128,N); %xb = DCT matrix
Hxb = zeros(N,N);
for i = 1:N
    Hx       = conv(xb(:,i),hrf);
    Hxb(:,i) = Hx(k + 128);
end
xb = xb(129:end,:);

% Get confounds (in scan time) and constant term
%--------------------------------------------------------------------------
X0 = xY(1).X0;
M  = size(X0,2);

% Specify covariance components; assume neuronal response is white
% treating confounds as fixed effects
%--------------------------------------------------------------------------
Q = speye(N,N)*N/trace(Hxb'*Hxb);
Q = blkdiag(Q, speye(M,M)*1e6  );


% Get whitening matrix (NB: confounds have already been whitened)
%--------------------------------------------------------------------------
W = SPM.xX.W(Sess.row,Sess.row);


% Create structure for spm_PEB
%--------------------------------------------------------------------------
clear P
P{1}.X = [W*Hxb X0];        % Design matrix for lowest level
P{1}.C = speye(N,N)/4;      % i.i.d assumptions
P{2}.X = sparse(N + M,1);   % Design matrix for parameters (0's)
P{2}.C = Q;

% Deconvolution
%=========================================================================
tic
% Iterate over voxels and deconvolve
total = 153594; counter = 0; percentage = 0;
x= SPM.xY.VY(1).dim(1);
y= SPM.xY.VY(1).dim(2);
z= SPM.xY.VY(1).dim(3);
t = SPM.nscan;
deconvData = zeros(x,y,z,t);

% cd into analyze_dir to be able to calculate beta
cd(analyze_dir);
for ix = 1:x
    for iy = 1:y
        for iz = 1:z
           
            thisVoxel = spm_data_read(VI,'xyz',[ix;iy;iz]);
            thisVoxel = spm_filter(SPM.xX.K,W*thisVoxel);
            
            % Remove null space of contrast
            beta  = spm_data_read(SPM.Vbeta,'xyz',[ix;iy;iz]);
            if ~isnan(mean(beta)) 
                thisVoxel = thisVoxel - spm_FcUtil('Y0',SPM.xCon(xY.Ic),SPM.xX.xKXs,beta);
            end
            % Simple deconvolution
            % -------------------------------------------------------------
            C  = spm_PEB(thisVoxel,P);
            xn = xb*C{2}.E(1:N);
            xn = spm_detrend(xn);
            
            % Save variables (NOTE: xn is in microtime and does not account for
            % slice timing shifts). To convert to BOLD signal convolve with a hrf.
            % Use a microtime to scan time index to convert to scan time: e.g.,
            % k = 1:NT:N*NT; where NT = number of bins per TR = TR/dt or SPM.xBF.T
            % and N = number of scans in the session. Finally account for slice
            % timing effects by shifting the index accordingly.
            %----------------------------------------------------------------------
            deconvVoxel = downsample(xn,16);  
            deconvData(ix, iy, iz,:) = deconvVoxel;
            
            counter = counter +1;
            percentage = (counter / total ) * 100;
           
        end
    end
end
toc

% CD back into working folder to save file
cd(pwd);
save([subj.curSubj '_deconvolvedData.mat'], 'deconvData');



end









