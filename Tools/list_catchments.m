function [Switches] = list_catchments(Switches)
%
% [Switches] = list_catchments(Switches)
% 
% Makes available the catchment names to HOOPLA via the Switch structure
% 
% Input : 
%   - Switches: information about specified computation
% Output: 
%   - Switches: information about specified computation
%
% Programmed by A. Thiboult (2017)


% Load catchments names
tmp=load(fullfile('.','Data',Switches.timeStep ,'Misc','catchment_names.mat'));
Switches.nameC=tmp.nameC;