
%% master script to input images from Tiff stack, segment, and extract data
%assumes tiff script are in order by channel (not by cycle)

%% inputs

%inputs that change with every batch
imageDirectory = '/n/scratch2/aet13/FI17/FI17_registered/';
saveDirectory = '/n/scratch2/aet13/FI17/FI17_output/';
experiment = 'FI17_';
mag = '';
maxCycle = 10; 
FOVlimits = [1,1,5,5,3,7,1,3]; % array of  timepoint, row, column, and field start/stop

timepoints = {"03h", "24h", "48h", "72h"};
rows ={"A", "B", "C", "D", "E", "F", "G", "H"};
columns ={"01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11"};

%inputs that change with every experiment
channelNames = { 'Hoechst1', 'Hoechst2', 'Hoechst3',...
    'Hoechst4','Hoechst5', 'Hoechst6', 'Hoechst7',...
    'Hoechst8','Hoechst9','Hoechst10',...
    'pH_GFP_Z', 'pH_GFP_slice', 'GFP_inact', 'p_p38', 'p_TBK1',...
    'Nos2', 'Arg-1', 'Lamp1','Ki67','none_488',...
    'mTB_mCherry_Z','mTB_mCherry_slice', 'p_STAT1', 'p_ERK', 'CD68',...
    'cPARP', 'HMGB1', 'actin', 'cCas3', 'none_555',...
    'cy5none', 'p50', 'LC3', 'CD80', 'RelA', 'p_mTOR', 'p_STAT3',...
    'PDL1', 'p_STAT6', 'IRF1'};

bugGFP = {'pH_GFP_Z', 'pH_GFP_slice'}; %channel number bug GFP signal
bugmCherry = {'mTB_mCherry_Z','mTB_mCherry_slice'}; %channel number bug mCherry signal
punctaChannels = {'p_TBK1', 'LC3'}; %channels to be analyzed for puncta

%%prepare for batch on O2 - create list of unique timepoint and well
elements = {timepoints(FOVlimits(1):FOVlimits(2)),...
    rows(FOVlimits(3):FOVlimits(4)),...
    columns(FOVlimits(5):FOVlimits(6))};
combinations = cell(1, numel(elements)); %set up the varargout result
[combinations{:}] = ndgrid(elements{:});
combinations = cellfun(@(x) x(:), combinations,'uniformoutput',false);
jobs = [combinations{:}];

parallel = true;
CycIF_execute(jobs, parallel, imageDirectory, saveDirectory,...
    experiment, mag, maxCycle, FOVlimits, channelNames,...
    bugGFP, bugmCherry, punctaChannels)




