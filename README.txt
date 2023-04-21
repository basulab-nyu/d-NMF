README

1. System Requirements
	Software: MATLAB
		Necessary Toolboxes
			Image Processing Toolbox v11.3
			Statistics and Machine Learning Toolbox v12.1
			Bioinformatics Toolbox v4.15.1
	Versions Tested On: Windows 11, MATLAB R2021a
	Required non-standard Hardware: None

2. Installation guide
	Instructions
		No formal installation required
		Clone or download the github repository to your local machine
		Add the directory and its subfolders to your MATLAB path
	Typical install time: <5 minutes

3. Demo
	Instructions to run on data
		Sample datasets are available at 
			https://doi.org/10.5281/zenodo.7851333
		Primary script is MAIN_Run_DNMF.m
			Change the "files" variable on line 4 to point to the tiff file to be 				operated on
				This tiff file should be a motion corrected tiff stack, ideally 					sampled at 1-2 Hz.
			Options may be modified in lines 7-32
			Run MAIN_Run_DNMF.m
		"Fitness trace" script is MAIN_fitness_detect_skeleton.m
			Change files on lines 2-3 to point to the tiff file and associated 					DNMF_Out.m file, respectively
	Expected output
		MAIN_Run_DNMF.m
			During runtime, a figure will be displayed showing the progress of 				segmentation
			2 files will be written to the folder containing the tiff stack
				DNMF_Out.mat: Contains the spatial footprints in the variable 					"cROIs," the temporal components in the variable "Cs," as well 					as various other statistics. cROIs stores the ROIs in column 					form. Individual ROIs can be viewed using the command 					imagesc(reshape(cROIs(:,i),[szA szB])), where i is the ROI 					number to view, and szA and szB define the pixel size of the 					field of view (typically 512 x 512)
				RoiSet_Auto.zip: Contains ROI boundaries in a format readable 				by ImageJ
			Expected run time: 5-10 minutes for a 512 x 512 movie with 1000 				frames
		MAIN_fitness_detect_skeleton.m
			Script will generate 2 cell array variables:
				START_STOP: Start and stop indices for calcium transients for 					each ROI
				PEAK: Indices of the peak value of calcium transients for each 					ROI
4. Instructions for use
	Similar as above