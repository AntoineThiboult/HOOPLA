function ini_calibration(Switches, iC, iM, iE, iS)
%
% function ini_calibration(Switches, iC, iM, iE)
%
% Initialize calibration
%
% Input:
%   - Switches: user specification for calibration
%   - iC: catchment number
%   - iM: model number
%   - iE: PET number
%   - iS: SAR number
%
% Programmed by A. Thiboult (2016)

%#ok<*NASGU>
if ~exist(fullfile('Results','Calibration',Switches.timeStep),'dir')
    mkdir(fullfile('Results','Calibration',Switches.timeStep));
end

if exist(fullfile('Results','Calibration',Switches.timeStep,...
        sprintf('C%s_H%s_E%s_S%s.mat',Switches.nameC{iC},Switches.nameM{iM,1},Switches.nameE{iE},Switches.nameS{iS})),'file')~=2 || ...
        Switches.overWrite.on == 1
    
    %% Processing input data
    if Switches.verb.on; dispstat(sprintf('Initializing the calibration for catchment %s and %s...',Switches.nameC{iC},Switches.nameM{iM,1}),'keepprev');end
    
    % Data specification for catchment / parameters
    DataPath.dataObs=fullfile('Data',Switches.timeStep,'Hydromet_obs',strcat('Hydromet_obs_',Switches.nameC{iC}));
    DataPath.modelParamBound=fullfile('Data',Switches.timeStep,'Model_parameters','model_param_boundaries');
    DataPath.snowModelParamBound=fullfile('Data',Switches.timeStep,'Model_parameters','snow_model_param_boundaries');
    DataPath.dataMetFcast=fullfile('Data',Switches.timeStep,'Det_met_fcast',strcat('Met_fcast_',Switches.nameC{iC}));
    
    % Check that all necessary data are provided
    [DataObs] = checkData(Switches, DataPath);
    
    % Crop observed data according to specified dates and warm up
    if Switches.warmUpCompute.on == 1      % Warm up
        [DataObs, ~ , DataWu] = cropData(Switches, DataObs, iM, iE, iS);
    elseif Switches.warmUpCompute.on == 0  % No Warm up
        [DataObs] =             cropData(Switches, DataObs, iM, iE, iS);
    end
    
    %% Launch calibration
    if Switches.snowmeltCompute.on == 1
        if Switches.warmUpCompute.on == 1
            [Result, DataCal, SarResult] = calibration(Switches, DataObs, DataPath, iC, iM, iE, iS, DataWu);
        elseif Switches.warmUpCompute.on==0
            [Result, DataCal, SarResult] = calibration(Switches, DataObs, DataPath, iC, iM, iE, iS);
        end
        if Switches.verb.on; dispstat(sprintf('Saving results...\n\n\n'),'keepprev');end
        Result=structfun(@single,Result,'UniformOutput', false);
        DataCal=structfun(@single,DataCal,'UniformOutput', false);
        SarResult=structfun(@single,SarResult,'UniformOutput', false);
        save(fullfile('Results','Calibration',Switches.timeStep,sprintf('C%s_H%s_E%s_S%s.mat',Switches.nameC{iC},Switches.nameM{iM,1},Switches.nameE{iE},Switches.nameS{iS})),'Result','DataCal','SarResult','Switches')
        
    elseif Switches.snowmeltCompute.on == 0
        
        if Switches.warmUpCompute.on==1
            [Result, DataCal] = calibration(Switches, DataObs, DataPath, iC, iM, iE, iS, DataWu);
        elseif Switches.warmUpCompute.on==0
            [Result, DataCal] = calibration(Switches, DataObs, DataPath, iC, iM, iE, iS);
        end
        if Switches.verb.on; dispstat(sprintf('Saving results...\n\n\n'),'keepprev');end
        Result2=structfun(@single,Result,'UniformOutput', false);
        DataCal=structfun(@single,DataCal,'UniformOutput', false);
        save(fullfile('Results','Calibration',Switches.timeStep,sprintf('C%s_H%s_E%s_S%s.mat',Switches.nameC{iC},Switches.nameM{iM,1},Switches.nameE{iE},Switches.nameS{iS})),'Result','DataCal','Switches')
    end
else
     dispstat(sprintf('Simulation for catchment %s, %s, %s, and %s already exist.\nConsider setting the "overwrite" switch to 1.',...
        Switches.nameC{iC},Switches.nameM{iM,1},Switches.nameE{iE},Switches.nameS{iS}),'keepprev');
end

