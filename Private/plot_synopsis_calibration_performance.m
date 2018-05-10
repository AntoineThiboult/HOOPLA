% script to plot:
%   - boxplot by model
%   - boxplot by catchment
%   - empirical cumulative distribution of performance
%
% The score that is displayed is the one chosen as the objective function

close all
clear
clc

%% Import paths
% Result path
[importPath] = uigetdir('..','Please select the folder that contains your calibration files');
if ~importPath; error('No directory specified. plot_synopsis_calibration aborted'); end

% Make deterministic scores available
addpath(fullfile('..','Tools','Calibration'));

%% Identify catchment and model names
% List calibration files
listFile=dir(fullfile(importPath,'*.mat'));
listFileName={listFile.name};
if isempty(listFile); error('The selected folder doesn''t contain any .mat file'); end
% List of hydrological models
listModName=regexp(listFileName,'(?<=H).*(?=\_E)','match');
listModName=unique([listModName{:}]);
[~, reIndex] = sort( str2double( regexp( listModName, '\d+', 'match', 'once' ))); % Natural sorting
listModName=listModName(reIndex);
% List of catchments
listCatchName=regexp(listFileName,'(?<=C).*(?=_H)','match');
listCatchName=unique([listCatchName{:}]);

% Preallocation
perfCal = NaN(numel(listCatchName),numel(listModName));

%% Retrieve data
for iC=1:numel(listCatchName)
    for iM=1:numel(listModName)
        fileToLoad = dir(fullfile(importPath, sprintf('*C%s_*H%s_*.mat',listCatchName{iC}, listModName{iM})));
        if ~isempty(fileToLoad)
            tmp = load(fullfile(importPath, fileToLoad.name));
            perfCal(iC,iM) = tmp.Result.bestf;
        end
    end
end

%% Boxplot by model
figure
boxplot(perfCal,{1:size(perfCal,2)})
grid on
title('Performance by model (over catchments)')
set(gca,'XTickLabel',sprintf('%s\n',listModName{:}),'XTickLabelRotation', 90);
xlabel('Models'); ylabel(tmp.Switches.calibration.score)

%% Boxplot by catchment
figure
boxplot(perfCal', {1:size(perfCal,1)})
grid on
title('Performance by catchment (over models)')
set(gca,'XTickLabel',sprintf('C %s\n',listCatchName{:}),'XTickLabelRotation', 90);
xlabel('Catchment');ylabel(tmp.Switches.calibration.score)

%% Plot cumulative density function of performance
figure
ecdf(perfCal(:))
grid on
title('Cumulative distribution of performance')
xlabel(tmp.Switches.calibration.score); ylabel('Cumulative probability')