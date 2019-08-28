% ELECTROME_DIR points to Electrome TMS-EEG Networks
setenv('ELECTROME_DIR',fullfile(getenv('SCRATCH'),'Electrome'))

% DATADIR points to CC TMS-fMRI preprocessed ROI timeseries
[~,uname_str] = unix('uname');
if(strfind(uname_str,'Darwin'))
    setenv('CC_DATADIR',['' ...
                    'tms-fMRI/CC/roitimeseries/']);
elseif(strfind(getenv('SHERLOCK'),'2'))
    setenv('CC_DATADIR', ...
        fullfile(getenv('PI_SCRATCH'), ...
                'COMET',...
                'CausalConnectome',...
                'PPI'));  
    setenv('CC_DATADIR_ALT', ...
        fullfile(getenv('PI_SCRATCH'), ...
                'COMET',...
                'CausalConnectome',...
                'derivatives','denoiser'));  
end

if(strfind(getenv('SHERLOCK'),'2'))
    NETSCIDIR = fullfile(getenv('HOME'),'MATLAB');
    addpath(fullfile(NETSCIDIR,'netsci','external','matlab-cliquer'));
else
    NETSCIDIR = fullfile(getenv('HOME'),'MATLAB','netsci');
    addpath(fullfile(NETSCIDIR,'external','matlab-cliquer'));
end
addpath(NETSCIDIR) 
addpath(fullfile(NETSCIDIR,'netsci'));
addpath(fullfile(NETSCIDIR,'netsci','netsci'));
%run(fullfile(NETSCIDIR,'netsci','setup.m'))
addpath(fullfile(NETSCIDIR,'netsci','plot'));
addpath(fullfile(NETSCIDIR,'netsci','utils'));

GGMDIR = fullfile(getenv('HOME'),'MATLAB','ggmClass');
addpath(GGMDIR)
addpath(fullfile(GGMDIR,'solvers','QUIC'))
run(fullfile(GGMDIR,'setup.m'));

unix('module load math')
unix('module load openblas')