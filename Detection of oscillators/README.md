# General info for ‘Detection of oscillators’

The code is largely based on the Gaussian Process approach developed by Dr Nick Phillips that is found at https://github.com/ManchesterBioinference/GPosc and includes routines developed by Nick which are integrated with the GPML toolbox by C Rasmussen and N Hannes (2010) "Gaussian processes for machine learning (GPML) toolbox" Journal of Machine Learning Research, 11:3011 3015. Please unarchive GPML.zip before running.

The routines were further customised by Dr Veronica Biga to include Hilbert reconstructions and detection of fold change by pairing peaks and subsequent troughs (Manning et al. (2019)) as well as a paired false discovery rate analysis where multiple embryos analysed in the same conditions can be considered, as detailed in Soto & Biga et al. (2020).

Upon using this code, please cite:

Biga V & Soto X (2020) “Dynamic properties of noise and Her6 levels are optimised by miR-9, allowing the decoding of the Her6 oscillator”, EMBO Journal 39:e103558.
Manning CS et al. (2019) “Quantitative single-cell live imaging links HES5 dynamics with cell-state and fate in murine neurogenesis”, Nature Communications (10), 2835.
Phillips NE et al (2017) "Identifying stochastic oscillations in single-cell live imaging time series using Gaussian processes", PLoS Comput Biol 13(5): e1005479.

Readme for detection of oscillators
There are two main files, for which the use is detailed below:
1)	Her6_main_exportHilbert.m allows running the analysis of oscillators for each track. To select a particular dataset from our data, uncomment the content of one section for example %%HV Venus norm found in lines 8-13.
Notes: To run on new datasets, you will need to:
-organise the data into an .xls or .csv file with time as the first vector followed by track1, track 2 etc; you will also need to collect 1-3  ‘background tracks’ by tracking an area of the dish with no cells;
-set ‘fnames’ in line 9 to your datafile name;
-set ‘exptnames’ in line 10 to a simulation name similar to fnames;
-set ‘coldata’ to the start and end columns corresponding to tracks that need to be analysed; it is recommend to exclude tracks that are too short; do not include any null columns;
- set ‘bkgrdata’ to the columns containing background tracks;
-set ‘strow’ to the row corresponding to time 0; when reading files with headers, Matlab may include information that is a header rather than part of a track-this then excludes any spurious values not part of the track.
  
Outcome- when the run is complete the code should generate a folder with the name as defined in ‘exptnames’ and this should contain various images and files. The important files area the ‘xx nofdr.mat’ file containing all the simulation variables and the ‘SummaryFile.xls’ where the statistics are summarised. The stats are final with the exception of the False Discovery Rate tests, see 2).  
2)	getFDRbyexpt.m carries out a paired false discovery rate test with 3% statistical significance. It is designed to consider a pair of HV and HVP embryos imaged and analysed in the same way, representing 1 experimental repeat.
  Notes: To run this on new data, load the new simulation files generated in 1) in lines 6 and 12 correspondingly. The significance level of the FDR can be adjusted in line 87.
Outcome- when complete, this routine will display the % oscillators based on a 3% FDR threshold and generate a globalFDR xx.mat file (the name for this is defined and can be adjusted in line 100). It will also export a 0 or 1 variable in column H of each ‘SummaryFile.xls’: 1 indicates oscillatory activity and 0 indicates aperiodic.


