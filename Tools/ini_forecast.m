function ini_forecast(Switches, iC, iM, iE, iS)
%
% function ini_forecast(Switches, iC, iM, iE)
%
% Initialize forecast
%
% Input:
%   - Switches: user specification for forecast
%   - iC: catchment number
%   - iM: model number
%   - iE: PET number
%   - iS: SAR number
%
% Programmed by A. Thiboult (2016)

%#ok<*NASGU>
%#ok<*STRNU>

if ~exist(fullfile('Results','Forecast',Switches.timeStep),'dir')
    mkdir(fullfile('Results','Forecast',Switches.timeStep));
end

if exist(fullfile('Results','Forecast',Switches.timeStep,...
        sprintf('C%s_H%s_E%s_S%s.mat',Switches.nameC{iC},Switches.nameM{iM,1},Switches.nameE{iE},Switches.nameS{iS})),'file')~=2 || ...
        Switches.overWrite.on == 1
    
    if Switches.forecast.metEns.on == 0 % Deterministic weather forecast
        
        %% Processing input data
        if Switches.verb.on; dispstat(sprintf('Initializing the forecast for catchment %s and %s...',Switches.nameC{iC},Switches.nameM{iM,1}),'keepprev');end
        
        % Data specification for catchment / parameters
        DataPath.dataObs=fullfile('Data',Switches.timeStep,'Hydromet_obs',strcat('Hydromet_obs_',Switches.nameC{iC}));
        DataPath.dataMetFcast=fullfile('Data',Switches.timeStep,'Det_met_fcast',strcat('Met_fcast_',Switches.nameC{iC}));
        DataPath.modelParam=fullfile('Data',Switches.timeStep,'Model_parameters','Calibrated_param');
        
        % Check that all necessary data are provided
        [DataObs, DataMetFcast] = checkData(Switches, DataPath);
        
        % Crop observed data according to specified dates and warm up
        if Switches.warmUpCompute.on == 1      % Warm up
            [DataObs, DataMetFcast , DataWu] = cropData(Switches, DataObs, iM, iE, iS, DataMetFcast);
        elseif Switches.warmUpCompute.on == 0  % No Warm up
            [DataObs, DataMetFcast] =          cropData(Switches, DataObs, iM, iE, iS, DataMetFcast);
        end
        
        %% Launch forecast
        if Switches.snowmeltCompute.on == 1
            if Switches.warmUpCompute.on == 1
                [Result, DataFcast, SarResult] = forecast(Switches, DataObs, DataMetFcast, DataPath, iC, iM, iE, iS, DataWu);
            elseif Switches.warmUpCompute.on==0
                [Result, DataFcast, SarResult] = forecast(Switches, DataObs, DataMetFcast , DataPath, iC, iM, iE, iS);
            end
            if Switches.verb.on; dispstat(sprintf('Saving results...\n\n\n'),'keepprev');end
            Result=structfun(@single,Result,'UniformOutput', false); Result.DateFcast=double(Result.DateFcast);
            DataFcast=structfun(@single,DataFcast,'UniformOutput', false);
            SarResult=structfun(@single,SarResult,'UniformOutput', false);
            save(fullfile('Results','Forecast',Switches.timeStep,sprintf('C%s_H%s_E%s_S%s.mat',Switches.nameC{iC},Switches.nameM{iM,1},Switches.nameE{iE},Switches.nameS{iS})),'Result','DataFcast','SarResult','Switches')
            
        elseif Switches.snowmeltCompute.on == 0
            
            if Switches.warmUpCompute.on==1
                [Result, DataFcast] = forecast(Switches, DataObs, DataMetFcast, DataPath, iC, iM, iE, iS, DataWu);
            elseif Switches.warmUpCompute.on==0
                [Result, DataFcast] = forecast(Switches, DataObs, DataMetFcast, DataPath, iC, iM, iE, iS);
            end
            if Switches.verb.on; dispstat(sprintf('Saving results...\n\n\n'),'keepprev');end
            Result=structfun(@single,Result,'UniformOutput', false); Result.DateFcast=double(Result.DateFcast);
            DataFcast=structfun(@single,DataFcast,'UniformOutput', false);
            save(fullfile('Results','Forecast',Switches.timeStep,sprintf('C%s_H%s_E%s_S%s.mat',Switches.nameC{iC},Switches.nameM{iM,1},Switches.nameE{iE},Switches.nameS{iS})),'Result','DataFcast','Switches')
        end
        
        
    elseif Switches.forecast.metEns.on == 1 % Meteorological ensemble prediction system
        
        %% Initialization of Result matrices for MEPS
        % Data Pre-check
        %         if Switches.forecast.perfectFcast.on == 1
        %             error('%s%s','Perfect forecast and ensemble meteorological forcing are incompatible. ',...
        %                 'It would be highly inefficient to run several times the exact same meteorological forecast. ',...
        %                 'Consider disabling either the perfect forecast switch or the meteorological ensemble forecast switch')
        %         end
        metFile = matfile(fullfile('Data',Switches.timeStep,'Ens_met_fcast',sprintf('Met_fcast_%s.mat',Switches.nameC{iC})));
        nbMetMb       = numel(metFile.Met_fcast);
        if nbMetMb == 0
            error('No meteorological ensemble forecast for catchment %s', Switches.nameC{iC});
        end
        
        % Preallocation
        nbTimeStep    = numel(datenum(Switches.fcastStart) : str2double( regexprep(Switches.timeStep,'h','') )/24 : datenum(Switches.fcastEnd));
        if Switches.DA.on == 0; Switches.DA.N = 1; end
        Result.Qfcast = NaN(nbMetMb, nbTimeStep, Switches.forecast.hor, Switches.DA.N,'single');
        Result.Qs     = NaN(nbMetMb, nbTimeStep, Switches.DA.N,'single');
        if Switches.snowmeltCompute.on == 1
            SarResult.runoffD  = NaN(nbMetMb, nbTimeStep, Switches.DA.N,'single');
        end
        
        for iMetMb=1:nbMetMb
            
            %% Processing input data
            if Switches.verb.on; dispstat(sprintf('Initializing the forecast for catchment %s, %s, and meteorological member %i...',Switches.nameC{iC},Switches.nameM{iM,1}, iMetMb),'keepprev');end
            
            % Data specification for catchment / parameters
            DataPath.dataObs=fullfile('Data',Switches.timeStep,'Hydromet_obs',strcat('Hydromet_obs_',Switches.nameC{iC}));
            DataPath.dataMetFcast=fullfile('Data',Switches.timeStep,'Ens_met_fcast',sprintf('Met_fcast_%s',Switches.nameC{iC}));
            DataPath.dataMetFcastEnsMb=iMetMb;
            DataPath.modelParam=fullfile('Data',Switches.timeStep,'Model_parameters','Calibrated_param');
            
            % Check that all necessary data are provided
            [DataObs, DataMetFcast] = checkData(Switches, DataPath);
            
            % Crop observed data according to specified dates and warm up
            if Switches.warmUpCompute.on == 1      % Warm up
                [DataObs, DataMetFcast , DataWu] = cropData(Switches, DataObs, iM, iE, iS, DataMetFcast);
            elseif Switches.warmUpCompute.on == 0  % No Warm up
                [DataObs, DataMetFcast] =          cropData(Switches, DataObs, iM, iE, iS, DataMetFcast);
            end
            
            %% Launch forecast
            if Switches.snowmeltCompute.on == 1
                if Switches.warmUpCompute.on == 1
                    [iMetMbResult, iMetMbDataFcast, iMetMbSarResult] = forecast(Switches, DataObs, DataMetFcast, DataPath, iC, iM, iE, iS, DataWu);
                elseif Switches.warmUpCompute.on==0
                    [iMetMbResult, iMetMbDataFcast, iMetMbSarResult] = forecast(Switches, DataObs, DataMetFcast , DataPath, iC, iM, iE, iS);
                end
                Result.Qfcast(iMetMb,:,:,:)   = iMetMbResult.Qfcast;
                Result.Qs(iMetMb,:,:)         = iMetMbResult.Qs;
                SarResult.runoffD(iMetMb,:,:)  = iMetMbSarResult.runoffD;
                
            elseif Switches.snowmeltCompute.on == 0
                if Switches.warmUpCompute.on==1
                    [iMetMbResult, iMetMbDataFcast] = forecast(Switches, DataObs, DataMetFcast, DataPath, iC, iM, iE, iS, DataWu);
                elseif Switches.warmUpCompute.on==0
                    [iMetMbResult, iMetMbDataFcast] = forecast(Switches, DataObs, DataMetFcast, DataPath, iC, iM, iE, iS);
                end
                Result.Qfcast(iMetMb,:,:,:)   = iMetMbResult.Qfcast;
                Result.Qs(iMetMb,:,:)         = iMetMbResult.Qs;
            end
        end
        
        %% Save results
        DataFcast=iMetMbDataFcast;
        DataFcast=structfun(@single,DataFcast,'UniformOutput', false);
        Result.DateFcast=iMetMbResult.DateFcast;
        if Switches.verb.on; dispstat(sprintf('Saving results...\n\n\n'),'keepprev');end
        if Switches.snowmeltCompute.on == 1
            save(fullfile('Results','Forecast',Switches.timeStep,sprintf('C%s_H%s_E%s_S%s.mat',Switches.nameC{iC},Switches.nameM{iM,1},Switches.nameE{iE},Switches.nameS{iS})),'Result','DataFcast','SarResult','Switches','-v7.3')
        elseif Switches.snowmeltCompute.on == 0
            save(fullfile('Results','Forecast',Switches.timeStep,sprintf('C%s_H%s_E%s_S%s.mat',Switches.nameC{iC},Switches.nameM{iM,1},Switches.nameE{iE},Switches.nameS{iS})),'Result','DataFcast','Switches','-v7.3')
        end
        
    end
else
    dispstat(sprintf('Simulation for catchment %s, %s, %s, and %s already exist.\nConsider setting the "overwrite" switch to 1.',...
        Switches.nameC{iC},Switches.nameM{iM,1},Switches.nameE{iE},Switches.nameS{iS}),'keepprev');
end
