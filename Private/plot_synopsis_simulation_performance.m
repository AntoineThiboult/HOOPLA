% script to plot:
%   - boxplot by model
%   - boxplot by catchment
%   - empirical cumulative distribution of performance
%
% Please specify the idMetric variable to choose the score. 1=RMSE, 2=RSMEsqrt, etc.

close all
clear
clc

%% Choice of the metric [EDIT idMetric]
% idMetric:         1        2        3         4      5        6       7     8       9       10       11     12       13     14       15       16       17
metricName=     {'RMSE' 'RMSEsqrt' 'RMSElog' 'MSE' 'MSEsqrt' 'MSElog' 'MAE' 'NSE' 'NSEsqrt' 'NSEinv' 'PVE' 'PVEabs' 'Balance' 'r'     'bKGE'  'gKGE'   'KGEm'} ;
idMetric = 8;


%% Import paths
% Result path
[importPath] = uigetdir('..','Please select the folder that contains your simulation files');
if ~importPath; error('No directory specified. plot_synopsis_simulation aborted'); end

% Make deterministic scores available
addpath(fullfile('..','Tools','Calibration'));

%% Identify catchment and model names
% List Simulation files
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
perfSim = NaN(numel(listCatchName),numel(listModName));

%% Retrieve data
for iC=1:numel(listCatchName)
    for iM=1:numel(listModName)
        fileToLoad = dir(fullfile(importPath, sprintf('*C%s_*H%s_*.mat',listCatchName{iC}, listModName{iM})));
        if ~isempty(fileToLoad)
            tmp = load(fullfile(importPath, fileToLoad.name));
            score = det_scores(tmp.DataSim.Q, mean(tmp.Result.Qs, 2));
            perfSim(iC,iM) = score(idMetric);
        end
    end
end

%% Boxplot by model
figure
boxplot(perfSim,{1:size(perfSim,2)})
grid on
title('Performance by model (over catchments)')
set(gca,'XTickLabel',sprintf('%s\n',listModName{:}),'XTickLabelRotation', 90);
xlabel('Models'); ylabel(metricName{idMetric})

%% Boxplot by catchment
figure
boxplot(perfSim', {1:size(perfSim,1)})
grid on
title('Performance by catchment (over models)')
set(gca,'XTickLabel',sprintf('C %s\n',listCatchName{:}),'XTickLabelRotation', 90);
xlabel('Catchment');ylabel(metricName{idMetric})

%% Plot cumulative density function of performance
figure
ecdf(perfSim(:))
grid on
title('Cumulative distribution of performance')
xlabel(metricName{idMetric}); ylabel('Cumulative probability')