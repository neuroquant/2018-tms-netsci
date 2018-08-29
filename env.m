% ELECTROME_DIR points to Electrome TMS-EEG Networks
setenv('ELECTROME_DIR',fullfile(getenv('SCRATCH'),'Electrome'))

% DATADIR points to CC TMS-fMRI preprocessed ROI timeseries
setenv('DATADIR',['/Volumes/MACBOOKUSB/Datasets/' ...
                    'tms-fMRI/CC/roitimeseries/']);                


NETSCIDIR = fullfile(getenv('HOME'),'MATLAB','netsci');
addpath(NETSCIDIR) 
addpath(fullfile(NETSCIDIR,'netsci'));
addpath(fullfile(NETSCIDIR,'netsci','plot'));
addpath(fullfile(NETSCIDIR,'netsci','utils'));
addpath(fullfile(NETSCIDIR,'external','matlab-cliquer')); 

GGMDIR = fullfile(getenv('HOME'),'MATLAB','ggmClass');
addpath(GGMDIR)
