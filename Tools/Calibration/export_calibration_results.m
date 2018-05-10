function export_calibration_results(varargin)
%
% export_calibration_results(varargin)
%
% Collect the calibration results in ~/Results/Calibration to create the
% structure Bestxs used in Simulation and Forecast. Can be used either with
% input (Switches) or independently without any argument (but then requires
% some confirmation from user)
%
% Input
%   varargin = Switches (see launch_hoopla for Switches details)
% 
% Output 
%   none: the results are saved to the
%         ./Data/timeStep/Model_parameters/Calibrated_param.mat file. 
%
% Programmed by A. Thiboult (2017)

%% Verification of inputs and user's preference / safety not to erase previous files
if isempty(varargin)
    % Select folder that contains the data to export
    [importPath] = uigetdir;    
    % Retrieve time step
    listCal = dir(fullfile(importPath,'*.mat'));
    if isempty(listCal); error('No calibration results found in selected directory'); end
    tmp = load(fullfile(importPath,listCal(1).name));
    exportPath=fullfile('Data',tmp.Switches.timeStep,'Model_parameters','Calibrated_param.mat');
    if exist(exportPath,'file') == 2
        % Confirmation
        prompt = '\nAre you sure to overwrite the existing calibrated model parameters ? Y = yes\n';
        overWriteAns = input(prompt,'s');
        if strcmpi(overWriteAns,'Y') == 1
            runExportCal(exportPath,importPath)
        else
            warning('Calibration export aborted')
        end
    else
       runExportCal(exportPath,importPath)
    end
elseif ~isempty(varargin)
    Switches=varargin{1};
    exportPath=fullfile('Data',Switches.timeStep,'Model_parameters','Calibrated_param.mat');
    importPath=fullfile('Results','Calibration',Switches.timeStep);
    if exist(exportPath,'file') ~= 2 || Switches.overWrite.on == 1        
        runExportCal(exportPath,importPath)
    else
        warning('Calibration export aborted. Consider setting the switch Switches.overWrite.on=1')
    end
end
end

%% Perform calibration parameters exportation
function runExportCal(exportPath,importPath)
listCal=dir(fullfile(importPath,'*.mat'));
for iCal=1:numel(listCal)
    tmpStr =regexp(listCal(iCal).name,'[_.]','split');
    cName = tmpStr{1}(2:end); % Catchment name
    mName = tmpStr{2}(2:end); % Model name
    calFile=load(fullfile(importPath,listCal(iCal).name));
    Bestxs.(sprintf('Bestx_%s',cName)).(sprintf('%s',mName))=calFile.Result.bestParam; % Parameters
    Bestxs.(sprintf('Bestx_%s',cName)).(sprintf('%s_CalSwitches',mName))=calFile.Switches; % Specificities of the calibration
end
save(exportPath,'Bestxs')
fprintf('Done !\n')
end