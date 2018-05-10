function [calScore]=objFunction(modelParam,userdata)
%
% function [calScore]=objFunction(modelParam,userdata)
%
% Objective function that is called by the calibration algorithm. Perform
% simulation with the tested parameter set and provide score to minimize.
% All scores are transformed to ensure that they are meant to be minimized.
% 
% Input: 
%   modelParam: model parameters that are evaluated
%   userdata: all hydrological data necessary for calibration
% 
% Output: 
%   calScore: objective function value
% 
% Programmed by A. Thiboult (2016)

%% Warm Up
if userdata.Switches.warmUpCompute.on == 1
    % Launch warm up
    if userdata.Switches.snowmeltCompute.on == 1
        [ParamWu, SarParamWu] = warmUp(userdata.Switches, userdata.DataWu, userdata.iniHydroModel,...
            userdata.hydroModel, userdata.iniPetModel, userdata.petModel, userdata.iniSarModel, userdata.sarModel, modelParam);
    else
        [ParamWu]              = warmUp(userdata.Switches, userdata.DataWu, userdata.iniHydroModel,...
            userdata.hydroModel, userdata.iniPetModel, userdata.petModel, userdata.iniSarModel, userdata.sarModel, modelParam);
    end
end

%% Model initialization for calibration

% Compute potential evapotranspiration
if userdata.Switches.petCompute.on == 1
    [PetData] = userdata.iniPetModel(userdata.Switches, userdata.DataCal);
    [userdata.DataCal.E] = userdata.petModel(PetData);
end

% Snow accounting model initialization
if userdata.Switches.snowmeltCompute.on == 1
    [SarResult, SarParam] = userdata.iniSarModel(userdata.Switches, userdata.DataCal, modelParam);
end

% Hydrological model initialization
[Result,Param] = userdata.iniHydroModel(userdata.Switches, userdata.DataCal.Date, modelParam);

% Initialization of states with WarmUp
if userdata.Switches.warmUpCompute.on == 1
    Param(:) = ParamWu;
    if userdata.Switches.snowmeltCompute.on == 1
        SarParam(:) = SarParamWu;
    end
end

%% Run simulation
if userdata.Switches.snowmeltCompute.on == 1
    %% With snow accounting
    for t = 1 : length(userdata.DataCal.Date)
        [SarResult.runoffD(t), SarParam] = ...
            userdata.sarModel(userdata.DataCal.Pt(t), userdata.DataCal.T(t), userdata.DataCal.Tmax(t),...
            userdata.DataCal.Tmin(t), userdata.DataCal.Date(t,:), SarParam);
        [Result.Qs(t,1),Param] = userdata.hydroModel(SarResult.runoffD(t), userdata.DataCal.E(t), Param);
    end    
elseif userdata.Switches.snowmeltCompute.on == 0
    %% No Snow accounting
    for t = 1 : length(userdata.DataCal.Date)
        [Result.Qs(t,1),Param] = userdata.hydroModel(userdata.DataCal.Pt(t), userdata.DataCal.E(t), Param);
    end
end

%% Evaluation of the calibration
if userdata.Switches.calibration.rmWinter.on == 1
    idWinter=userdata.DataCal.Date(:,2)== 12 | userdata.DataCal.Date(:,2)== 1 | userdata.DataCal.Date(:,2)== 2 | userdata.DataCal.Date(:,2)== 3;
    [~, allCalScores]=det_scores(userdata.DataCal.Q(~idWinter), Result.Qs(~idWinter));
    calScore=allCalScores(userdata.idScore);
else
    [~, allCalScores]=det_scores(userdata.DataCal.Q, Result.Qs);
    calScore=allCalScores(userdata.idScore);
end