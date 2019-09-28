% This function initialises the subject variables for the batch analysis of
% the Movie Watching data.
% _________________________________________________________________________
% Last checked: September 2019

% Lorena Freitas
% $Id: initialize_vars.m 11 2016-12-01 16:26:24F Lorena $


function [b] = initialize_vars(subject,group)

% SPM info
%--------------------------------------------------------------------------
b.spmDir = fileparts(which('spm')); % path to SPM installation


% Directory information
%--------------------------------------------------------------------------
dataDir = '<GENERAL PATH FOR ALL SUBJECTS DATA>';


% Subject information
%--------------------------------------------------------------------------
b.group = group;
b.curSubj = subject; % subject code (e.g., 'sujet01')
b.dataDir = strcat(dataDir, filesep , b.curSubj, filesep); % make data directory subject-specific

% Data type to use (deconvolved = 1 on non-deconvolved = 0)
%--------------------------------------------------------------------------
b.deconvolveData = 1;

% Folder for functional files
%--------------------------------------------------------------------------
b.dir_fonc = '<PATH TO FUNCTIONAL DATA, e.g., /warped/ >';


end
