# README

This is an Etkinlab dataset acquired between 2013 to 2018. 
Description: Cross sectional dataset of structural and functional MRI in non-trauma exposed healthy, trauma exposed healthy, traume exposed patient populations. fMRI was acquired under multiple conditions including resting state, and task acquired on day 1 and interleaved TMS/fMRI to upto 17 cortical targets on day 2. 

func/ses-d2/
------------
Note that the first 3 volumes of the tms tasks are non-steady state. Volumes 4 to 167 correspond to the actual duration of the TMS stimulation task. 

### Comments from BIDS Curators ###
===================================
Jan 2019: Manjari Narayan
===================================

Currently a subset of fMRI conditions for the healthy population has been organized into BIDS format. 

Known Issues
 - Defacing has not been performed on structurals. Should be performed for anatomical T1w scans. 
 - `T1w.json` was created using dcm2niix whereas `_bold.json` was manually created by reading variables from ge header files using rdgehdr on the Pxxxx.7.hdr files. The information is assumed to be identical across all subjects and only differs across scan sequences used between day1 and day2 sessions. There are optional json fields that are not populated given the lack of pfile to nifti bids sidecar tools. 
 - Gary Glover's pfile reconstruction of spiral fMRI data into Pxxx.7.nii (no physiological denoising) and Pxxx.7.den.nii (physiological denoising) should eventually be provided in a sourcedata/ directory as described in the BIDS documentation.
 - Origin of current `_bold.nii.gz` files is unknown. In future iterations these files ought to have an additional qualifier recon-<recon-type> to explain which of of Gary Glover's Pfile to Pxxxx.7.nii reconstruction algorithms were used. Ideally we should have nifti images of both kinds.  
 - physiological recordings, extensive particiapant phenotype, behavioral, scan.tsv information are available for this dataset but have not yet been curated AFAIK. 
 - Behavioral data such as `_events.tsv` for fMRI tasks have not been made. Since the tms triggering is identical for all tasks, I have provided these. 
 - All the fMRI data with tms stimuli are named `task-singlepulse<BrainRegion>`. The corresponding json sidecars or event files should also contain targeting information such as the exact coordinates of TMS targetting and perhaps the spatial extent of the TMS electrical field.
