#* Author: Serafeim Loukas, EPFL, serafeim.loukas@epfl.ch, 17 July 2021
#* Script to overlay an activation map onto the spm152 structural volume, with 
#* specific predefined parameters and settings.

import gl, os

#* Where to save results
save_to = '/Users/loukas/Desktop/' # do not forget the last slash
path_to_caps = '/Users/loukas/Desktop/PPI_CAPs/Server_results/60%/caps'

#* Get file names
cap_files = [f for f in sorted(os.listdir(path_to_caps)) if f.startswith('PPI')&f.endswith('.nii') ]

#* Loop and plot
for f in range(len(cap_files)):
	path = path_to_caps + '/' + cap_files[f]

	gl.resetdefaults()
	
	#* Set Multi panel view mode
	# Display Axial (1), Coronal (2), Sagittal (4), Flipped Sagittal (8), MPR (16), Mosaic (32) or Rendering (64)
	#gl.view(32)
	#gl.mosaic("A L+ H -0.1 V -0.1 -14 2 18; 34 50 66 S X R");
	gl.view(16)
	
	#* Changes the background color, for example backcolor(255, 0, 0) will set a bright red background
	gl.backcolor(255, 255, 255)
	
	#* Open background image
	gl.loadimage('/Users/loukas/Desktop/PPI_CAPs/template/spm152_wrapped2neonatal.nii')
	gl.minmax(0, 40, 80)
	
	#* Open overlay: show positive regions
	gl.overlayloadsmooth(1) # no smoothing
	gl.overlayload(path)
	gl.colorname(1,"4hot")
	gl.minmax(1, 0.05, 2)
	#* Make the layer (0 for background, 1 for 1st overlay) transparent(0), translucent (~50) or opaque (100).
	gl.opacity(1,80)
	
	#open overlay: show negative regions
	gl.overlayloadsmooth(1)
	gl.overlayload(path)
	gl.minmax(2, -0.05, -2)
	#gl.colorname (2,"7cool")
	gl.colorname (2,"5winter")
	#* Make the layer (0 for background, 1 for 1st overlay) transparent(0), translucent (~50) or opaque (100).
	gl.opacity(2,80)
	
	#* Set colorbar position (0=off, 1=top, 2=right).
	gl.colorbarposition(1)
	#gl.orthoviewmm(17, 1, 25)
	gl.orthoviewmm(5, 1, 25)
	gl.savebmp(save_to + 'image_{}.png'.format(cap_files[f].split('_')[-3]))	


