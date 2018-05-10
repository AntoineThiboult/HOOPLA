% script to plot:
%   - boxplot by model
%   - boxplot by catchment
%   - empirical cumulative distribution of performance
%
% 1) Please specify the the lead time for which scores will be computed and
%    ploted. It should be specified in time steps. ex leadTime = [1 3 10];
% 2) Please specify the idMetric variable to choose the score. 1=RMSE, 2=RSMEsqrt, etc.

close all
clear
clc

%% User preferences [EDIT idMetric and listLeadTimes]

% Choice of the lead time to represent
listLeadTimes = [1 3 5];

% Choice of the metric
% idMetric:         1        2        3         4      5        6       7     8       9       10       11     12       13     14       15       16       17
metricName=     {'RMSE' 'RMSEsqrt' 'RMSElog' 'MSE' 'MSEsqrt' 'MSElog' 'MAE' 'NSE' 'NSEsqrt' 'NSEinv' 'PVE' 'PVEabs' 'Balance' 'r'     'bKGE'  'gKGE'   'KGEm'} ;
idMetric = 8;

%% Import paths
% Result path
[importPath] = uigetdir('..','Please select the folder that contains your forecast files');
if ~importPath; error('No directory specified. plot_synopsis_forecast aborted'); end

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
perfFcast = NaN(numel(listCatchName),numel(listModName),numel(listLeadTimes));

%% Retrieve data
for iC=1:numel(listCatchName)
    for iM=1:numel(listModName)
        fileToLoad = dir(fullfile(importPath, sprintf('*C%s_*H%s_*.mat',listCatchName{iC}, listModName{iM})));
        if ~isempty(fileToLoad)
            tmp = load(fullfile(importPath, fileToLoad.name));
            
            QfcastLT = mean(tmp.Result.Qfcast,4);
            if tmp.Switches.forecast.metEns.on == 1
                QfcastLT = squeeze(mean(QfcastLT,1));
            end
            
            for iLT=1:numel(listLeadTimes)
                QrefLT = cat(1,tmp.DataFcast.Q(listLeadTimes(iLT):end),NaN(listLeadTimes(iLT),1)); % Shifting the obs to match the forecast matrix (see HOOPLA manual, Sec. 6.3)
                score = det_scores(QrefLT, QfcastLT(:,listLeadTimes(iLT)));
                perfFcast(iC,iM,iLT) = score(idMetric);
            end
        end
    end
end

%% Boxplot by model
figure
for iLT=1:numel(listLeadTimes)
    subplot(1,numel(listLeadTimes),iLT)
    boxplot(squeeze(perfFcast(:,:,iLT)),{1:size(perfFcast,2)})
    grid on
    title(sprintf('Performance by model (over catchments)\n for lead time %i',listLeadTimes(iLT)))
    set(gca,'YLim',[min(perfFcast(:))-0.01, max(perfFcast(:))+0.01])
    set(gca,'XTickLabel',sprintf('%s\n',listModName{:}),'XTickLabelRotation', 90);
    xlabel('Models'); ylabel(metricName{idMetric})
end

%% Boxplot by catchment
figure
for iLT=1:numel(listLeadTimes)
    subplot(1,numel(listLeadTimes),iLT)
    boxplot(squeeze(perfFcast(:,:,iLT))',{1:size(perfFcast,1)})
    grid on
    title(sprintf('Performance by catchment (over models)\n for lead time %i',listLeadTimes(iLT)))
    set(gca,'YLim',[min(perfFcast(:))-0.01, max(perfFcast(:))+0.01])
    set(gca,'XTickLabel',sprintf('C %s\n',listCatchName{:}),'XTickLabelRotation', 90);
    xlabel('Catchment');ylabel(metricName{idMetric})
end

%% Plot cumulative density function of performance
figure
for iLT=1:numel(listLeadTimes)
    subplot(1,numel(listLeadTimes),iLT)
    ecdf(perfFcast(:,:,iLT))
    grid on
    title(sprintf('Cumulative distribution of performance\n for lead time %i',listLeadTimes(iLT)))
    set(gca,'XLim',[min(perfFcast(:))-0.01, max(perfFcast(:))+0.01])
    xlabel(metricName{idMetric}); ylabel('Cumulative probability')
end