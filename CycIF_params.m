%% master script to input images from Tiff stack, segment, and extract data
%assumes tiff script are in order by channel (not by cycle)

%% inputs

%inputs that change with every batch
imageDirectory = 'F:\FI11\Registered2\';
saveDirectory = 'C:\Users\Amy Thurber\Dropbox (Partners HealthCare)\Experiments\FI11\Matlab analysis\';
experiment = 'FI11_';
timepoint = '3h_';
mag = '20x_';
maxCycle = 5; 
r_start = 2; %first row to analyze
r_stop = 4; %last row to analyze
c_start = 2;%first column to analyze
c_stop = 7;%last column to analyze
f_start = 1;%first field to analyze
f_stop = 3; %last field to analyze

%inputs that change with every experiment
channelNames = { 'Hoechst1', 'Hoechst2', 'Hoechst3',...
    'Hoechst4','Hoechst5', 'pH_GFP',...
    'p_TBK1', 'p_cJun', 'Arg-1', 'p_p38', 'mTB_mCherry',...
    'p_ERK', 'p_JNK', 'CD68', 'p_STAT1',...
    'mitoRed', 'p_STAT6', 'p_mTOR', 'RelA', 'p_STAT3'};

bug_GFP = maxCycle + 1; %channel number bug GFP signal
bug_mCherry = (2*maxCycle) + 1; %channel number bug mCherry signal
puncta_channels = [3,4]; %channels to be analyzed for puncta
