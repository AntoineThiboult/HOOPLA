function [Result, DataSim, varargout] = forecast(Switches, DataSim, DataMetFcast, DataPath, iC, iM, iE, iS, varargin)
%
% function [Result, DataSim, SarResult] = forecast(Switches, DataSim, DataPath, iC, iM, iE, DataWu)
%
% Inputs:
%
%   Switches.
%       PET_comput      = 0 for provided PET (E in Data) ; 1 for external computation
%       Snowmelt_comput = 0 for no snowmelt ; 1 for needed snowmelt computation
%       Warm_up         = 0 for no warm-up ; 1 for needed warm-up (5 times mean Precip year)
%       ... (see launch_Hoopla for all Switches)
%   DataSim.
%       Date    = Simulated dates (yyyy/mm/dd/hh:mm:ss)
%       Q       = Observed streamfow (matrix size: nDay x 1)%
%       T       = Mean temperature (°C)
%       ... (see Doc/Hoopla_Manual for all fields)
%   DataMetFcast.
%       Date = date matrix for meteorological forecast (yyyy/mm/dd/hh:mm:ss)
%       Pt   = meteorological precipitation forecast
%       T    = meteorological temperature forecast
%       Tmin = meteorological precipitation forecast
%       Tmax = meteorological precipitation forecast
%       leadTime = vector that contains the forecast lead time in day unit
%   DataPath.
%       dataObsPath     = Data file name (ex : 'Data.mat') path
%       modelParamBound = model parameter boundaries path
%       snowModelParamBound = snow model parameter boundaries path
%   iC     = catchment number
%   iM     = model number
%   iE     = potential evapotranspiration formula number
%   DataWu = Data for warm up
%
%
% Outputs:
%
%   Results.
%       Qfcast    = Streamflow forecast with a size of
%                   - (nb of fcast issues, horizon, number of DA member) with data assimilation
%                   - (nb of fcast issues, horizon) in open loop fashion
%       Qs        = Simulated streamflow
%       DateFcast = Date matrix
%   DataSim.
%       Date    = Simulated dates (yyyy/mm/dd/hh:mm:ss)
%       Q       = Observed streamfow (matrix size: nDay x 1)
%       E       = Potential evapotranspiration (mm). Note that the potential
%                 evapotranspiration may be computed but doesn't appear in
%                 the Result structure. It is treated as a data.
%       T       = Mean temperature (°C)
%       ... (see Doc/Hoopla_Manual for all fields)
%       Note that the suffix RP can be added to any value to denote random
%       perturbations
%   SarResult.
%       runoffD = runoff including snowmelt
%       Pg       = sarModel's solid precipitations (mm)
%       Pl       = sarModel's liquid precipitations (mm)
%       G        = sarModel's snowpack (mm)
%       snowMelt = sarModel's snowmelt (°C)
%
% Programmed by A. Thiboult (2016)

%#ok<*RHSFN> Code analyser doesn't recognize handles

%% Load calibrated parameters
load(DataPath.modelParam)
try
    modelParam=Bestxs.(sprintf('Bestx_%s',Switches.nameC{iC})).(Switches.nameM{iM,1});
catch
    error('No calibrated parameters available for catchment %s and model %s. Consider performing calibration.',...
        Switches.nameC{iC},Switches.nameM{iM,1})
end

%% Function handles
iniHydroModel=str2func(strcat('ini_',Switches.nameM{iM,1}));  % handle of the function ini_model iM
hydroModel=str2func(Switches.nameM{iM,1});                    % handle of the function model iM
if Switches.petCompute.on == 1
    iniPetModel=str2func(strcat('ini_',Switches.nameE{iE}));    % handle of the function ini_PET iE
    petModel=str2func(Switches.nameE{iE});                      % handle of the function PET iE
else
    iniPetModel = [];
    petModel = [];
end
if Switches.snowmeltCompute.on == 1
    iniSarModel=str2func(strcat('ini_',Switches.nameS{iS}));    % handle of the function ini_PET iE
    sarModel=str2func(Switches.nameS{iS});                      % handle of the function PET iE
else
    iniSarModel = [];
    sarModel = [];
end
if Switches.DA.on == 1
    daModel=str2func(strcat(Switches.DA.tech,'_func'));         % handle of the assimilation function
end

%% Reservoirs to update
if Switches.DA.on == 1
    allModelUpdatedRes=load(fullfile('Data',Switches.timeStep,'Misc','reservoir_to_update'));
    Switches.DA.updatedRes=allModelUpdatedRes.(Switches.nameM{iM,1});
end

%% Warm Up
if Switches.warmUpCompute.on == 1
    % Input argument
    DataWu=varargin{1};
    % Launch warm up
    if Switches.snowmeltCompute.on == 1
        [ParamWu, SarParamWu] = warmUp(Switches, DataWu, iniHydroModel, hydroModel, iniPetModel, petModel, iniSarModel, sarModel, modelParam);
    else
        [ParamWu]              = warmUp(Switches, DataWu, iniHydroModel, hydroModel, iniPetModel, petModel, iniSarModel, sarModel, modelParam);
    end
end

%% Data assimilation initialization
if Switches.DA.on == 1
    [DataSim, w]=ini_DA(Switches,DataSim);% Noisy inputs
end

%% Simulation
if Switches.verb.on && ~Switches.parallelCompute.on; dispstat('Beginning of the simulation...','keepprev'); dispstat('Simulation progress: 0%');end

%% Simulation with Data Assimilation
if Switches.DA.on == 1
    %% Initialization of models for the Data Assimilation
    
    % Compute potential evapotranspiration
    if Switches.petCompute.on == 1
        DataSim.ERP=NaN(size(DataSim.Date,1),Switches.DA.N);
        selectField  ={'TpetRP','TminRP','TmaxRP'};
        selectField2 ={'T'     ,'Tmin'  ,'Tmax'  };
        for imDA=1:Switches.DA.N
            for iField=1:numel(selectField)
                DataSimPETtmp.(selectField2{iField})=DataSim.(selectField{iField})(:,imDA);
            end
            DataSimPETtmp = mergeStruct(DataSim, DataSimPETtmp); % Retrieve data from DataSim and merge them with DataDA
            [PetData] = iniPetModel(Switches,DataSimPETtmp);
            [DataSim.ERP(:,imDA)] = petModel(PetData); % Compute PET for the ith DA member
        end
    end
    
    % Snow accounting model initialization
    if Switches.snowmeltCompute.on == 1
        [SarResult, SarParam] = iniSarModel(Switches, DataSim, modelParam);
        SarResult = structfun(@(x) ( permute(x, [1,3,2]) ), SarResult, 'UniformOutput', false); % add a ghost dimension
        SarResult = structfun(@(x) ( repmat(x, [1, Switches.DA.N, 1]) ), SarResult, 'UniformOutput', false); % replicate for the matrix for the N DA members
        SarParam(1:Switches.DA.N)=SarParam; % replicate for the structure for the N DA members
    end
    
    % Hydrological model initialization
    [Result,Param] = iniHydroModel(Switches, DataSim.Date, modelParam);
    Result = structfun(@(x) ( permute(x, [1,3,2]) ), Result, 'UniformOutput', false); % add a ghost dimension
    Result = structfun(@(x) ( repmat(x, [1, Switches.DA.N, 1]) ), Result, 'UniformOutput', false); % replicate for the parameters for the N DA members
    Param(1:Switches.DA.N)=Param; % replicate for the structure for the N DA members
    
    % Initialization of states with WarmUp
    if Switches.warmUpCompute.on == 1
        Param(:) = ParamWu;
        if Switches.snowmeltCompute.on == 1
            SarParam(:) = SarParamWu;
        end
    end
    
    % Initialization forecast matrices
    nbFcastIssue            = length(DataMetFcast.Date);
    DataMetFcast.E          = NaN(nbFcastIssue,Switches.forecast.hor);
    Result.Qfcast           = NaN(nbFcastIssue,Switches.forecast.hor,Switches.DA.N);
    if Switches.snowmeltCompute.on == 1
        SarResultFcast.runoffD  = NaN(nbFcastIssue,Switches.forecast.hor,Switches.DA.N);
    end
    
    %% Run simulation
    if Switches.snowmeltCompute.on == 1
        %% With snow accounting
        for t = 1 : size(DataSim.Date,1)
            
            if Switches.verb.on && ~mod(t,floor(size(DataSim.Date,1)/20)) && ~Switches.parallelCompute.on; dispstat(sprintf('Simulation progress: %d %%',round(100*t/size(DataSim.Date,1)))); end
            
            for imDA=1:Switches.DA.N
                [SarResult.runoffD(t,imDA), SarParam(imDA)] = ...
                    sarModel(DataSim.PtRP(t,imDA), DataSim.TsnoRP(t,imDA),...
                    DataSim.TmaxRP(t,imDA), DataSim.TminRP(t,imDA), DataSim.Date(t,:), SarParam(imDA));
                [Result.Qs(t,imDA), Param(imDA)] = ...
                    hydroModel(SarResult.runoffD(t,imDA), DataSim.ERP(t,imDA), Param(imDA));
            end
            
            % Perform DA
            if ~mod(t,Switches.DA.dt)
                if ~any(isnan(DataSim.QRP(t,:)))
                    switch Switches.DA.tech
                        case 'EnKF'
                            [Param]   = daModel(Param, Result.Qs(t,:), DataSim.Q(t,:), DataSim.QRP(t,:), DataSim.eQRP(t,:), Switches.DA, w);
                        case 'PF'
                            [Param, w]= daModel(Param, Result.Qs(t,:), DataSim.Q(t,:), DataSim.QRP(t,:), DataSim.eQRP(t,:), Switches.DA, w);
                        case 'PertOnly' % Do nothing.
                    end
                end
            end
            
            % Hydrological forecast
            if ~any(isnan(DataMetFcast.Pt(t,:)))
                if Switches.petCompute.on == 1 % PET
                    DataMetFcastPET.Date = datevec(datenum(DataMetFcast.Date(t,:))+DataMetFcast.leadTime);
                    DataMetFcastPET.T    = DataMetFcast.T(t,:)';
                    DataMetFcastPET.Tmin = DataMetFcast.Tmin(t,:)';
                    DataMetFcastPET.Tmax = DataMetFcast.Tmax(t,:)';
                    DataMetFcastPET      = mergeStruct(DataSim, DataMetFcastPET); % Retrieve data from DataSim and merge them with DataFcastPET
                    [PetData]           = iniPetModel(Switches,DataMetFcastPET);
                    [DataMetFcast.E(t,:)]= petModel(PetData); % Compute PET for the ith DA member
                end
                % Loop over DA members
                for imDA=1:Switches.DA.N
                    % Set the initial forecast states with the ones obtained from simulation
                    SarParamFcast=SarParam;
                    ParamFcast=Param;
                    % Loop over lead times
                    for iLt=1:Switches.forecast.hor
                        % Snow
                        iDate=datevec(datenum(DataMetFcast.Date(t,:))+DataMetFcast.leadTime(iLt));
                        [SarResultFcast.runoffD(t,iLt,imDA), SarParamFcast(imDA)] = ...
                            sarModel(DataMetFcast.Pt(t,iLt), DataMetFcast.T(t,iLt),...
                            DataMetFcast.Tmax(t,iLt), DataMetFcast.Tmin(t,iLt), iDate, SarParamFcast(imDA));
                        % Hydrological model
                        [Result.Qfcast(t,iLt,imDA), ParamFcast(imDA)] = ...
                            hydroModel(SarResultFcast.runoffD(t,iLt,imDA), DataMetFcast.E(t,iLt), ParamFcast(imDA));
                    end
                end
            end
        end
        
        % Output
        varargout{1}=SarResult;
        
    elseif Switches.snowmeltCompute.on == 0
        %% No Snow accounting
        for t = 1 : size(DataSim.Date,1)
            
            if Switches.verb.on && ~mod(t,floor(size(DataSim.Date,1)/20)) && ~Switches.parallelCompute.on; dispstat(sprintf('Simulation progress: %d %%',round(100*t/size(DataSim.Date,1)))); end
            
            for imDA=1:Switches.DA.N
                [Result.Qs(t,imDA), Param(imDA)] = ...
                    hydroModel(DataSim.PtRP(t,imDA), DataSim.ERP(t,imDA), Param(imDA));
            end
            
            % Perform DA
            if ~mod(t,Switches.DA.dt)
                if ~any(isnan(DataSim.QRP(t,:)))
                    switch Switches.DA.tech
                        case 'EnKF'
                            [Param]   = daModel(Param, Result.Qs(t,:), DataSim.Q(t,:), DataSim.QRP(t,:), DataSim.eQRP(t,:), Switches.DA, w);
                        case 'PF'
                            [Param, w]= daModel(Param, Result.Qs(t,:), DataSim.Q(t,:), DataSim.QRP(t,:), DataSim.eQRP(t,:), Switches.DA, w);
                        case 'PertOnly' % Do nothing.
                    end
                end
            end
            
            % Hydrological forecast
            if ~any(isnan(DataMetFcast.Pt(t,:)))
                if Switches.petCompute.on == 1 % PET
                    DataMetFcastPET.Date = datevec(datenum(DataMetFcast.Date(t,:))+DataMetFcast.leadTime);
                    DataMetFcastPET.T    = DataMetFcast.T(t,:)';
                    DataMetFcastPET.Tmin = DataMetFcast.Tmin(t,:)';
                    DataMetFcastPET.Tmax = DataMetFcast.Tmax(t,:)';
                    DataMetFcastPET      = mergeStruct(DataSim, DataMetFcastPET); % Retrieve data from DataSim and merge them with DataFcastPET
                    [PetData]           = iniPetModel(Switches,DataMetFcastPET);
                    [DataMetFcast.E(t,:)]= petModel(PetData); % Compute PET for the ith DA member
                end
                % Loop over DA members
                for imDA=1:Switches.DA.N
                    % Set the initial forecast states with the ones obtained from simulation
                    ParamFcast=Param;
                    % Loop over lead times
                    for iLt=1:Switches.forecast.hor
                        % Hydrological model
                        [Result.Qfcast(t,iLt,imDA), ParamFcast(imDA)] = ...
                            hydroModel(DataMetFcast.Pt(t,iLt), DataMetFcast.E(t,iLt), ParamFcast(imDA));
                    end
                end
            end
        end
    end
end


%% Open Loop simulation
if Switches.DA.on == 0
    %% Initialization of models for the open loop
    
    % Compute potential evapotranspiration
    if Switches.petCompute.on == 1
        [PetData] = iniPetModel(Switches,DataSim);
        [DataSim.E] = petModel(PetData); % Compute PET for the ith DA member
    end
    
    % Snow accounting model initialization
    if Switches.snowmeltCompute.on == 1
        [SarResult, SarParam] = iniSarModel(Switches, DataSim, modelParam);
    end
    
    % Hydrological model initialization
    [Result,Param] = iniHydroModel(Switches, DataSim.Date, modelParam);
    
    % Initialization of states with WarmUp
    if Switches.warmUpCompute.on == 1
        Param(:) = ParamWu;
        if Switches.snowmeltCompute.on == 1
            SarParam(:) = SarParamWu;
        end
    end
    
    % Initialization matrices forecast
    nbFcastIssue            = length(DataMetFcast.Date);
    DataMetFcast.E          = NaN(nbFcastIssue,Switches.forecast.hor);
    Result.Qfcast           = NaN(nbFcastIssue,Switches.forecast.hor);
    if Switches.snowmeltCompute.on == 1
        SarResultFcast.runoffD  = NaN(nbFcastIssue,Switches.forecast.hor);
    end
    
    %% Run simulation
    if Switches.snowmeltCompute.on == 1
        %% With snow accounting
        for t = 1 : size(DataSim.Date,1)
            
            if Switches.verb.on && ~mod(t,floor(size(DataSim.Date,1)/20)) && ~Switches.parallelCompute.on; dispstat(sprintf('Simulation progress: %d %%',round(100*t/size(DataSim.Date,1)))); end
            
            [SarResult.runoffD(t), SarParam] = ...
                sarModel(DataSim.Pt(t), DataSim.T(t), DataSim.Tmax(t),...
                DataSim.Tmin(t), DataSim.Date(t,:), SarParam);
            [Result.Qs(t,1),Param] = hydroModel(SarResult.runoffD(t), DataSim.E(t), Param);
            
            % Hydrological forecast
            if ~any(isnan(DataMetFcast.Pt(t,:)))
                if Switches.petCompute.on == 1 % PET
                    DataMetFcastPET.Date = datevec(datenum(DataMetFcast.Date(t,:))+DataMetFcast.leadTime);
                    DataMetFcastPET.T    = DataMetFcast.T(t,:)';
                    DataMetFcastPET.Tmin = DataMetFcast.Tmin(t,:)';
                    DataMetFcastPET.Tmax = DataMetFcast.Tmax(t,:)';
                    DataMetFcastPET      = mergeStruct(DataSim, DataMetFcastPET); % Retrieve data from DataSim and merge them with DataFcastPET
                    [PetData]           = iniPetModel(Switches,DataMetFcastPET);
                    [DataMetFcast.E(t,:)]= petModel(PetData); % Compute PET for day t and all lead times
                end
                % Set the initial forecast states with the ones obtained from simulation
                SarParamFcast=SarParam;
                ParamFcast=Param;
                % Loop over lead times
                for iLt=1:Switches.forecast.hor
                    % Snow
                    iDate=datevec(datenum(DataMetFcast.Date(t,:))+DataMetFcast.leadTime(iLt));
                    [SarResultFcast.runoffD(t,iLt), SarParamFcast] = ...
                        sarModel(DataMetFcast.Pt(t,iLt), DataMetFcast.T(t,iLt),...
                        DataMetFcast.Tmax(t,iLt), DataMetFcast.Tmin(t,iLt), iDate, SarParamFcast);
                    % Hydrological model
                    [Result.Qfcast(t,iLt), ParamFcast] = ...
                        hydroModel(SarResultFcast.runoffD(t,iLt), DataMetFcast.E(t,iLt), ParamFcast);
                end
            end
        end
        
        % Output
        varargout{1}=SarResult;
        
    elseif Switches.snowmeltCompute.on == 0
        %% No Snow accounting
        for t = 1 : size(DataSim.Date,1)
            
            if Switches.verb.on && ~mod(t,floor(size(DataSim.Date,1)/20)) && ~Switches.parallelCompute.on; dispstat(sprintf('Simulation progress: %d %%',round(100*t/size(DataSim.Date,1)))); end
            
            [Result.Qs(t,1),Param] = hydroModel(DataSim.Pt(t), DataSim.E(t), Param);
            
            % Hydrological forecast
            if ~any(isnan(DataMetFcast.Pt(t,:)))
                if Switches.petCompute.on == 1 % PET
                    DataMetFcastPET.Date = datevec(datenum(DataMetFcast.Date(t,:))+DataMetFcast.leadTime);
                    DataMetFcastPET.T    = DataMetFcast.T(t,:)';
                    DataMetFcastPET.Tmin = DataMetFcast.Tmin(t,:)';
                    DataMetFcastPET.Tmax = DataMetFcast.Tmax(t,:)';
                    DataMetFcastPET      = mergeStruct(DataSim, DataMetFcastPET); % Retrieve data from DataSim and merge them with DataFcastPET
                    [PetData]           = iniPetModel(Switches,DataMetFcastPET);
                    [DataMetFcast.E(t,:)]= petModel(PetData); % Compute PET for day t and all lead times
                end
                % Set the initial forecast states with the ones obtained from simulation
                ParamFcast=Param;
                % Loop over lead times
                for iLt=1:Switches.forecast.hor
                    % Hydrological model
                    [Result.Qfcast(t,iLt), ParamFcast] = ...
                        hydroModel(DataMetFcast.Pt(t,iLt), DataMetFcast.E(t,iLt), ParamFcast);
                end
            end
        end
    end
end
Result.DateFcast = repmat(datenum(DataMetFcast.Date),[1, numel(DataMetFcast.leadTime)]) + repmat(DataMetFcast.leadTime,[size(DataMetFcast.Date,1),1]);