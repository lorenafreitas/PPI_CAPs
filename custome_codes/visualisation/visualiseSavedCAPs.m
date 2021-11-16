%
% Modified by Serafeim Loukas, July 19, 2021
%

%ppicapsDir = '/Users/loukas/Desktop/PPI_CAPs/Server_results/60%/';
%save_caps = [ppicapsDir 'caps/'];

ppicapsDir = '/Users/loukas/Desktop/PPI_CAPs/results/';
ppicapsDir = '/Users/loukas/Desktop/';
save_caps = ppicapsDir;

initialPath = save_caps;
savePath = initialPath;
caps =  spm_select('FPList', initialPath, '^PPI-CAPs_2Level_K5.*\.nii$');
if isempty(caps); error("No CAPs were found") ;end

frame = 'axial'; % coronal / sagittal / axial
for cap = 1:size(caps,1)
    img     = job_Visualise_hot_winter(caps(cap,:), frame);
    [filepath,name,ext] = fileparts(caps(cap,1:(end-4)));
    print(img, [savePath name '_' frame],'-dpdf','-r300','-fillpage');
    close all;
end
