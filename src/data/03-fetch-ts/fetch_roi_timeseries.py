'''
    File name: fetch_roi_timeseries.py
    Author: Manjari Narayan
    Date created: 01/15/2019
    Date last modified: 10/31/2017
    Python Version: 3.6
    Description: Script to extract roi level timeseries after fmriprep and denoiser
    Project: Causal Connectome, 2018-tms-netsci
'''

from __future__ import division

from nilearn.input_data import NiftiMasker, NiftiMapsMasker, NiftiLabelsMasker
from nilearn.connectome import ConnectivityMeasure
from matplotlib import pyplot as plt
from scipy.signal import lfilter
from operator import itemgetter
from collections import Counter
from datetime import datetime
from itertools import groupby
from nilearn import datasets
import nilearn.signal
import nibabel as nib
import nilearn.image
import pandas as pd
import numpy as np
import argparse
import nilearn
import sys
import errno
import os


subject = os.environ.get('SUBJECT')
print(subject)
#subject = 'NTHC1001'
cleandir = os.path.join(os.environ.get('PI_SCRATCH'),
                        "COMET",
                        "CausalConnectome",
                        "derivatives/fmriprep-fsl/denoiser/",
                        "sub-%s"%subject)

keys = os.listdir(cleandir)
#keys = np.array([x.split("_bold") for x in keys if 'task' in x]).flatten()
keys = np.unique([x for x in keys if 'task' and '.nii.gz' in x]).tolist()

# prepdir = os.environ.get('PREPDIR')
# subprep = os.path.join(prepdir,"sub-"+subject,"MNINonLinear/Results")

print(datetime.now().strftime("%a %b %d %H:%M:%S"))
print("Creating roi timeseries")


for key in keys:
    print("extracting session "+key)

    #prepfile = os.path.join(subprep,key,key+".nii.gz") #original file
    #totaltp = nib.load(prepfile).shape[3]
    # if totaltp <= 10:
    #     continue

    # imgfile = os.path.join(cleandir,key,key+'_removed_first10_despiked_masked_mvmreg%s_cmpc_bp.nii.gz'%gsr)
    imgfile = os.path.join(cleandir,key)
    
    #####################################
    #  Gordon, Schaefer, Buckner atlas #
    #####################################

    atlasfile = os.path.join(os.environ.get("PI_HOME"),"resources",
             'SchaeferYeo2018/MNI/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz')
    atlasfile2 = os.path.join(os.environ.get("PI_HOME"),"resources",
              'Gordon333/Parcels_MNI_222.nii')
    atlasfile3 = os.path.join(os.environ.get("PI_HOME"),"resources",
              'SchaeferYeo2018/MNI/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz')
    atlasfile4 = os.path.join(os.environ.get("PI_HOME"),"resources","parcellations","MNI","2mm",
              'Shen268.nii.gz')          
    atlasfile5 = os.path.join(os.environ.get("PI_HOME"),"resources",
    'SchaeferYeo2018/MNI/Schaefer2018_300Parcels_7Networks_order_FSLMNI152_2mm.nii.gz')
    subcort_atlasfile = os.path.join(os.environ.get("PI_HOME"),"resources",
 'Choi_JNeurophysiol12_MNI152/Choi2012_7Networks_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz')
    cerebellum_atlasfile = os.path.join(os.environ.get("PI_HOME"),"resources",
 'Buckner_JNeurophysiol11_MNI152/Buckner2011_7Networks_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz')

    # extract signals
    masker = NiftiLabelsMasker(labels_img=atlasfile,smoothing_fwhm=4,standardize=False,detrend=False,low_pass=None,high_pass=None,verbose=5)
    masker2 = NiftiLabelsMasker(labels_img=atlasfile2,smoothing_fwhm=4,standardize=False,detrend=False,low_pass=None,high_pass=None,verbose=5)
    masker3 = NiftiLabelsMasker(labels_img=atlasfile3,smoothing_fwhm=4,standardize=False,detrend=False,low_pass=None,high_pass=None,verbose=5) 
    masker4 = NiftiLabelsMasker(labels_img=atlasfile4,smoothing_fwhm=4,standardize=False,detrend=False,low_pass=None,high_pass=None,verbose=5)        
    masker5 = NiftiLabelsMasker(labels_img=subcort_atlasfile,smoothing_fwhm=4,standardize=False,detrend=False,low_pass=None,high_pass=None,verbose=5)
    masker6 = NiftiLabelsMasker(labels_img=cerebellum_atlasfile,smoothing_fwhm=4,standardize=False,detrend=False,low_pass=None,high_pass=None,verbose=5)
    masker7 = NiftiLabelsMasker(labels_img=atlasfile5,smoothing_fwhm=4,standardize=False,detrend=False,low_pass=None,high_pass=None,verbose=5)

    # save parcellated time series
    try:
        os.mkdir(os.path.join(cleandir,key.replace('.nii.gz','')+'_roits'))
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass
        
    ####### Atlas 1 ########
    time_series = masker.fit_transform(imgfile, 
                                        confounds=None)
    outfile = os.path.join(cleandir,key.replace('.nii.gz','')+'_roits',
                                    key.replace('.nii.gz','')+"_Schaefer200_Yeo7Networks.csv")
    np.savetxt(outfile,time_series)
    ####### Atlas 2 ########
    time_series2 = masker2.fit_transform(imgfile,
                                        confounds=None)
    outfile2 = os.path.join(cleandir,key.replace('.nii.gz','')+'_roits',
                                    key.replace('.nii.gz','')+"_Gordon333.csv")
    np.savetxt(outfile2,time_series2)
    ####### Atlas 3 ########
    time_series3 = masker3.fit_transform(imgfile,
                                        confounds=None)
    outfile3 = os.path.join(cleandir,key.replace('.nii.gz','')+'_roits',
                                    key.replace('.nii.gz','')+"_Schaefer100_Yeo7Networks.csv")
    np.savetxt(outfile3,time_series3)
    ####### Atlas 7 ########
    time_series7 = masker7.fit_transform(imgfile,
                                        confounds=None)
    outfile7 = os.path.join(cleandir,key.replace('.nii.gz','')+'_roits',
                                    key.replace('.nii.gz','')+"_Schaefer300_Yeo7Networks.csv")
    np.savetxt(outfile7,time_series7)
    
    ####### Atlas 4 ########
    time_series4 = masker4.fit_transform(imgfile,
                                        confounds=None)
    outfile4 = os.path.join(cleandir,key.replace('.nii.gz','')+'_roits',
                                    key.replace('.nii.gz','')+"_Shen268.csv")
    np.savetxt(outfile4,time_series4)
    
    ####### Atlas 5 ########
    time_series5 = masker5.fit_transform(imgfile,
                                        confounds=None)
    outfile5 = os.path.join(cleandir,key.replace('.nii.gz','')+'_roits',
                                    key.replace('.nii.gz','')+"_Choi7.csv")
    np.savetxt(outfile5,time_series5)
    
    ####### Atlas 6 ########
    time_series6 = masker6.fit_transform(imgfile,
                                        confounds=None)
    outfile6 = os.path.join(cleandir,key.replace('.nii.gz','')+'_roits',
                                    key.replace('.nii.gz','')+"_Buckner7.csv")
    np.savetxt(outfile6,time_series6)
    
    #####################################
    # Points to remember
    # smoothing_fwhm: float, optional
    #   If smoothing_fwhm is not None, it gives the full-width half maximum in millimeters of the spatial smoothing to apply to the signal.
    # 
    # resampling_target: {“data”, “labels”, None}, optional. 
    # The default option is that labels are resampled to the data. 

    # ####### Atlas Sub-cortical ########
    # try:
    #     time_series_subcort = subcortmasker.fit_transform(imgfile)
    #     outfile4 = os.path.join(cleandir,key.replace('.nii.gz','')+'_roits',
    #                                 key.replace('.nii.gz','')+"_Striatum_Yeo7Networks.csv")
    #     np.savetxt(outfile4,time_series_subcort)
    # except OSError as exc:
    #     print(exc.errno)
    #     print('Skipping Choi subcortical atlas')
    # ####### Atlas Cerebellum ########
    # try:
    #     time_series_cerebellum = cerebellummasker.fit_transform(imgfile)
    #     outfile5 = os.path.join(cleandir,key.replace('.nii.gz','')+'_roits',
    #                                 key.replace('.nii.gz','')+"_Cerebellum_Yeo7Networks.csv")
    #     np.savetxt(outfile5,time_series_cerebellum)
    # except OSError as exc:
    #     print(exc.errno)
    #     print('Skipping Cerebellar atlas')