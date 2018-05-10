function [Switches] = list_models(Switches)
%
% [Switches] = list_models(Switches)
%
% Makes available the models names (hydrological models, potential
% evapotranspiration formulae, and snow accounting routine)to HOOPLA via
% the Switch structure
%
% Input :
%   - Switches: information about specified computation
% Output:
%   - Switches: information about specified computation
%
% Programmed by A. Thiboult (2017)

% Load hydrological model names
tmp=load(fullfile('.','Data',Switches.timeStep ,'Misc','hydro_model_names.mat'));
Switches.nameM=tmp.nameM;

% Load potential evapotranspiration model names
if Switches.petCompute.on==1
    tmp=load(fullfile('.','Data',Switches.timeStep ,'Misc','pet_model_names.mat'));
    Switches.nameE=tmp.nameE;
else
    Switches.nameE={'NoPet'};
end

% Load snow accounting model names
if Switches.snowmeltCompute.on==1
    tmp=load(fullfile('.','Data',Switches.timeStep ,'Misc','snow_model_names.mat'));
    Switches.nameS=tmp.nameS;
else
    Switches.nameS={'NoSar'};
end