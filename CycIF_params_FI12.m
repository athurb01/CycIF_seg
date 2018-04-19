%% master script to input images from Tiff stack, segment, and extract data
%assumes tiff script are in order by channel (not by cycle)

%% inputs

%inputs that change with every batch
imageDirectory = 'C:\Users\Amy Thurber\Dropbox (Partners HealthCare)\Experiments\FI12\FI12_registered\';
saveDirectory = 'C:\Users\Amy Thurber\Dropbox (Partners HealthCare)\Experiments\FI12\matlab_output\';
experiment = 'FI12_';
timepoint = '24h_';
mag = '';
maxCycle = 5; 
FOVlimits = [2,2,2,2,1,2]; % array of row start/stop, column start/stop, field start/stop

%inputs that change with every experiment
channelNames = { 'Hoechst1', 'Hoechst2', 'Hoechst3',...
    'Hoechst4','Hoechst5',...
    'pH_GFP', 'pH_GFP1z', 'GFP_inact', 'fitcnone', 'p_TBK1',...
    'mTB_mCherry','mTB_mCherry1z', 'cy3_inact', 'CD68', 'STAT1',...
    'cy5none', 'p50', 'cy5_inact', 'LC3', 'RelA'};

bugGFP = maxCycle + 1; %channel number bug GFP signal
bugmCherry = (2*maxCycle) + 1; %channel number bug mCherry signal
punctaChannels = [10,19]; %channels to be analyzed for puncta

[morphOut, fluorescenceOut] = CycIF_head(imageDirectory, saveDirectory,...
    experiment, timepoint, mag, maxCycle, FOVlimits, channelNames,...
    bugGFP, bugmCherry, punctaChannels)