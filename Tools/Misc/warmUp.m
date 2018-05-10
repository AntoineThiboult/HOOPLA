function [Param, varargout]=warmUp(Switches, DataWu, iniHydroModel, hydroModel, iniPetModel, petModel, iniSarModel, sarModel, modelParam)
%
% [Param, SarParam]=warmUp(Switches, DataWu, iniHydroModel, hydroModel, petModel, modelParam)
%
% Provide coherent hydrological model state variable values regarding past
% conditions. These state variable values are used in calibration/simulation/forecast
% as initial conditions
% 
% Input: 
%
%   Switches.
%       PET_comput      = 0 for provided PET (E in Data) ; 1 for external computation
%       Snowmelt_comput = 0 for no snowmelt ; 1 for needed snowmelt computation
%       Warm_up         = 0 for no warm-up ; 1 for needed warm-up (5 times mean Precip year)
%       ... (see launch_Hoopla for all Switches)
%   DataWu.
%       Date    = Simulated dates (yyyy/mm/dd/hh:mm:ss)
%       Q       = Observed streamfow (matrix size: nDay x 1)%
%       T       = Mean temperature (°C)
%       ... (see Doc/Hoopla_Manual for all fields)
%   iniHydroModel: handle of function for hydrological model initialization
%   hydroModel: handle of the hydrological model function
%   iniPetModel: handle of function for PET model initialization
%   petModel: handle of the PET function
%   modelParam: calibrated parameters for the hydrological model (and snow model)
%
%
% Output:
%
%   Param: hydrological calibrated parameters and initial state varaible values
%   SarParam: snow accounting routine calibrated parameters and initial state varaible values
%
% Programmed by A. Thiboult (2016)

runMode=dbstack;
if Switches.verb.on && ~strcmp(runMode(2).name,'objFunction'); dispstat('Warm up computation...','keepprev');end

%% Hydrological model initialization
[Result,Param] = iniHydroModel(Switches, DataWu.Date, modelParam);

%% Snow accounting model initialization
if Switches.snowmeltCompute.on == 1
    [SarResult, SarParam] = iniSarModel(Switches, DataWu, modelParam);
end

%% PET
if Switches.petCompute.on == 1
    [PetData] = iniPetModel(Switches, DataWu);
    [DataWu.E] = petModel(PetData);
end

%% Simulation
if Switches.snowmeltCompute.on == 1     % Snow accounting
    if Switches.exportLight.on == 1
        % Run simulation
        for t = 1 : length(DataWu.Date)
            [SarResult.runoffD(t), SarParam] = sarModel(DataWu.Pt(t), DataWu.T(t), DataWu.Tmax(t), DataWu.Tmin(t), DataWu.Date(t,:), SarParam);
            [Result.Qs(t,1),Param] = hydroModel(SarResult.runoffD(t), DataWu.E(t), Param);
        end
    elseif Switches.exportLight.on == 0
        % Run simulation
        for t = 1 : length(DataWu.Date)
            [SarResult.runoffD(t), SarParam, SarResult.Pg(t,:), SarResult.Pl(t,:),...
                SarResult.G(t,:), SarResult.eTg(t,:), SarResult.snowMelt(t,:)] = ...
                sarModel(DataWu.Pt(t), DataWu.T(t), DataWu.Tmax(t), DataWu.Tmin(t), DataWu.Date(t,:), SarParam);
            [Result.Qs(t,1),Param, Result.inter(t,:), Result.interq(t,:)] = ...
                hydroModel(SarResult.runoffD(t), DataWu.E(t), Param);
        end
    end    
    % Outputs
    varargout{1}=SarParam;
    
elseif Switches.snowmeltCompute.on == 0 % No Snow accounting
    if Switches.exportLight.on == 1
        % Run simulation
        for t = 1 : length(DataWu.Date)
            [Result.Qs(t,1),Param] = hydroModel(DataWu.Pt(t), DataWu.E(t), Param);
        end
    elseif Switches.exportLight.on == 0
        % Run simulation
        for t = 1 : length(DataWu.Date)
            [Result.Qs(t,1),Param, Result.inter(t,:), Result.interq(t,:)] = ...
                hydroModel(DataWu.Pt(t), DataWu.E(t), Param);
        end
    end
end