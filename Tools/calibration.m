function [Result, DataCal, varargout] = calibration(Switches, DataCal, DataPath, ~, iM, iE, iS, varargin)
%
% function [Result, DataSim, SarResult] = calibration(Switches, DataCal, DataPath, iC, iM, iE, DataWu)
%
% Inputs:
%
%   Switches.
%       PET_comput      = 0 for provided PET (E in Data) ; 1 for external computation
%       Snowmelt_comput = 0 for no snowmelt ; 1 for needed snowmelt computation
%       Warm_up         = 0 for no warm-up ; 1 for needed warm-up (5 times mean Precip year)
%       ... (see launch_Hoopla for all Switches)
%   DataCal.
%       Date    = Simulated dates (yyyy/mm/dd/hh:mm:ss)
%       Q       = Observed streamfow (matrix size: nDay x 1)%
%       T       = Mean temperature (°C)
%       ... (see Doc/Hoopla_Manual for all fields)
%   DataPath.
%       dataObsPath     = Data file name (ex : 'Data.mat') path
%       modelParamBound = model parameter boundaries path
%       snowModelParamBound = snow model parameter boundaries path
%   iM     = model number
%   iE     = potential evapotranspiration formula number
%   DataWu = Data for warm up
%
%
% Outputs:
%
%   Results.
%       bestParam = identified model best parameters
%       bestf     = best cost function
%       allBestf  = best cost function for each iteration
%       Qs        = Simulated streamflow with best parameters
%       DateCal   = Date matrix
%       inter     = hydrological model internal values
%       interq    = hydrological model flow components
%   DataCal.
%       Date    = Simulated dates (yyyy/mm/dd/hh:mm:ss)
%       Q       = Observed streamfow (matrix size: nDay x 1)
%       E       = Potential evapotranspiration (mm). Note that the potential
%                 evapotranspiration may be computed but doesn't appear in
%                 the Result structure. It is treated as a data.
%       T       = Mean temperature (°C)
%       ... (see Doc/Hoopla_Manual for all fields)
%   SarResult.
%       runoffD = runoff including snowmelt
%       Pg       = sarModel's solid precipitations (mm)
%       Pl       = sarModel's liquid precipitations (mm)
%       G        = sarModel's snowpack (mm)
%       snowMelt = sarModel's snowmelt (°C)
%
% Programmed by A. Thiboult (2016)

%#ok<*RHSFN> Code analyser doesn't recognize handles

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

%% Calibration settings

% Parameter boundaries
modelParamBound=load(DataPath.modelParamBound);
snowModelParamBound=load(DataPath.snowModelParamBound);
if Switches.calibration.snowCal.on == 1
    S_min=[modelParamBound.(Switches.nameM{iM,1}).sMin, snowModelParamBound.(Switches.nameS{iS}).sMin];
    S_ini=[modelParamBound.(Switches.nameM{iM,1}).sIni, snowModelParamBound.(Switches.nameS{iS}).sIni];
    S_max=[modelParamBound.(Switches.nameM{iM,1}).sMax, snowModelParamBound.(Switches.nameS{iS}).sMax];
elseif Switches.calibration.snowCal.on == 0
    S_min=[modelParamBound.(Switches.nameM{iM,1}).sMin, snowModelParamBound.(Switches.nameS{iS}).default];
    S_ini=[modelParamBound.(Switches.nameM{iM,1}).sIni, snowModelParamBound.(Switches.nameS{iS}).default];
    S_max=[modelParamBound.(Switches.nameM{iM,1}).sMax, snowModelParamBound.(Switches.nameS{iS}).default];
end

% Scores for the objective function
score={'RMSE','RMSEsqrt','RMSElog','MSE','MSEsqrt','MSElog','MAE','NSE',...
    'NSEsqrt','NSEinv','PVE','PVEabs','Balance','r','bKGE','gKGE','KGEm'};
userdata.idScore= find(cellfun(@(s) ( strcmp(Switches.calibration.score, s)), score));
userdata.orientScore = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0 ,0, 1, 1, 1, 1];

% Transfer of model information to the optimizitation code
userdata.Switches=Switches;
userdata.iniHydroModel=iniHydroModel;
userdata.hydroModel=hydroModel;
if Switches.petCompute.on == 1
    userdata.iniPetModel=iniPetModel;
    userdata.petModel=petModel;
else 
    userdata.iniPetModel = [];
    userdata.petModel = [];
end
if Switches.snowmeltCompute.on == 1
    userdata.iniSarModel=iniSarModel;
    userdata.sarModel=sarModel;
else
    userdata.iniSarModel = [];
    userdata.sarModel = [];
end

userdata.verbose=Switches.verb.on;
userdata.DataCal=DataCal;
if nargin == 8; userdata.DataWu=varargin{1}; end

%% Calibration
if Switches.verb.on;    
    dispstat('Beginning of the calibration...','keepprev');
end

switch Switches.calibration.method
    case 'DDS'
        maxiter=Switches.calibration.maxiter;
        [bestParam,bestf,allBestf] = dds('objFunction', S_ini, S_min, S_max, maxiter, userdata);
    case 'SCE'
        ngs=Switches.calibration.SCE.ngs;
        [bestParam,bestf,allBestf] = sce_ua('objFunction', S_ini, S_min, S_max, ngs, userdata);
end


%% Run simulation with best parameters
% Warm Up
if Switches.warmUpCompute.on == 1
    % Launch warm up
    if Switches.snowmeltCompute.on == 1
        [ParamWu, SarParamWu] = warmUp(Switches, userdata.DataWu, iniHydroModel,...
            hydroModel, iniPetModel, petModel, iniSarModel, sarModel, bestParam);
    else
        [ParamWu]              = warmUp(Switches, userdata.DataWu, iniHydroModel,...
            hydroModel, iniPetModel, petModel, iniSarModel, sarModel, bestParam);
    end
end

% Compute potential evapotranspiration
if Switches.petCompute.on == 1
    [PetData]  = iniPetModel(Switches,DataCal);
    [DataCal.E] = petModel(PetData);
end

% Snow accounting model initialization
if Switches.snowmeltCompute.on == 1
    [SarResult, SarParam] = iniSarModel(Switches, DataCal, bestParam);
end

% Hydrological model initialization
[Result,Param] = iniHydroModel(Switches, DataCal.Date, bestParam);

% Initialization of states with WarmUp
if Switches.warmUpCompute.on == 1
    Param(:) = ParamWu;
    if Switches.snowmeltCompute.on == 1
        SarParam(:) = SarParamWu;
    end
end

% Run simulation
if Switches.snowmeltCompute.on == 1
    %% With snow accounting
    if  Switches.exportLight.on == 1
        for t = 1 : length(DataCal.Date)
            if Switches.verb.on && ~mod(t,floor(length(DataCal.Date)/20)) && ~Switches.parallelCompute.on; dispstat(sprintf('Simulation progress: %d %%',round(100*t/length(DataCal.Date)))); end
            [SarResult.runoffD(t), SarParam] = ...
                sarModel(DataCal.Pt(t), DataCal.T(t), DataCal.Tmax(t),...
                DataCal.Tmin(t), DataCal.Date(t,:), SarParam);
            [Result.Qs(t,1),Param] = hydroModel(SarResult.runoffD(t), DataCal.E(t), Param);
        end
    elseif Switches.exportLight.on == 0
        for t = 1 : length(DataCal.Date)
            if Switches.verb.on && ~mod(t,floor(length(DataCal.Date)/20)) && ~Switches.parallelCompute.on; dispstat(sprintf('Simulation progress: %d %%',round(100*t/length(DataCal.Date)))); end
            [SarResult.runoffD(t), SarParam, SarResult.Pg(t,:), SarResult.Pl(t,:),...
                SarResult.G(t,:), SarResult.eTg(t,:), SarResult.snowMelt(t,:)] = ...
                sarModel(DataCal.Pt(t), DataCal.T(t), DataCal.Tmax(t),...
                DataCal.Tmin(t), DataCal.Date(t,:), SarParam);
            [Result.Qs(t,1),Param, Result.inter(t,:),Result.interq(t,:)] = ...
                hydroModel(SarResult.runoffD(t), DataCal.E(t), Param);
        end
    end
    % Output
    varargout{1}=SarResult;
    
elseif Switches.snowmeltCompute.on == 0
    %% No Snow accounting
    if Switches.exportLight.on == 1
        for t = 1 : length(DataCal.Date)
            if Switches.verb.on && ~mod(t,floor(length(DataCal.Date)/20)) && ~Switches.parallelCompute.on; dispstat(sprintf('Simulation progress: %d %%',round(100*t/length(DataCal.Date)))); end
            [Result.Qs(t,1),Param] = hydroModel(DataCal.Pt(t), DataCal.E(t), Param);
        end
    elseif Switches.exportLight.on == 0
        for t = 1 : length(DataCal.Date)
            if Switches.verb.on && ~mod(t,floor(length(DataCal.Date)/20)) && ~Switches.parallelCompute.on; dispstat(sprintf('Simulation progress: %d %%',round(100*t/length(DataCal.Date)))); end
            [Result.Qs(t,1),Param, Result.inter(t,:),Result.interq(t,:)] = ...
                hydroModel(DataCal.Pt(t), DataCal.E(t), Param);
        end
    end
end

%% Outputs
Result.bestParam=bestParam;
Result.bestf=bestf;
Result.allBestf = allBestf;
Result.DateCal=DataCal.Date;