%% master script to input images from Tiff stack, segment, and extract data
%assumes tiff script are in order by channel (not by cycle)

%% inputs

%inputs that change with every batch
imageDirectory = 'C:\Users\Amy Thurber\Dropbox (Partners HealthCare)\Experiments\FI12\FI12_registered\';
saveDirectory = 'C:\Users\Amy Thurber\Dropbox (Partners HealthCare)\Experiments\FI12\matlab_output\';
experiment = 'FI12_';
mag = '';
maxCycle = 5; 
FOVlimits = [3,3,2,2,2,2,1,1]; % array of  timepoint, row, column, and field start/stop

timepoints = {"03h", "24h", "48h", "72h"};
rows ={"A", "B", "C", "D", "E", "F", "G", "H"};
columns ={"01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11"};

%inputs that change with every experiment
channelNames = { 'Hoechst1', 'Hoechst2', 'Hoechst3',...
    'Hoechst4','Hoechst5',...
    'pH_GFP', 'pH_GFP1z', 'GFP_inact', 'fitcnone', 'p_TBK1',...
    'mTB_mCherry','mTB_mCherry1z', 'cy3_inact', 'CD68', 'STAT1',...
    'cy5none', 'p50', 'cy5_inact', 'LC3', 'RelA'};

bugGFP = maxCycle + 1; %channel number bug GFP signal
bugmCherry = (2*maxCycle) + 1; %channel number bug mCherry signal
punctaChannels = [10,19]; %channels to be analyzed for puncta

%%create list of unique timepoint and well
elements = {timepoints(FOVlimits(1):FOVlimits(2)),...
    rows(FOVlimits(3):FOVlimits(4)),...
    columns(FOVlimits(5):FOVlimits(6))};
combinations = cell(1, numel(elements)); %set up the varargout result
[combinations{:}] = ndgrid(elements{:});
combinations = cellfun(@(x) x(:), combinations,'uniformoutput',false);
jobs = [combinations{:}];

%% if statement to call parcluster or regular analysis function
for j = 1:length(jobs)
    timepoint = jobs(j,1);
    row = jobs(j,2);
    column = jobs(j,3);

%else

    [morphOut, fluorescenceOut] = CycIF_head(imageDirectory, saveDirectory,...
        experiment, timepoint, row, column, mag, maxCycle, FOVlimits, channelNames,...
        bugGFP, bugmCherry, punctaChannels);
end
