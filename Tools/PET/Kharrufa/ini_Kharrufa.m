function [PetData] = ini_Kharrufa(~, Data)
% [PetData] = ini_Oudin(Data)
%
% Compute the prerequisite values that are necessary for the computation of Oudin
%
%   Inputs
%       Data.
%           Date   = Date / date (datenum)           
%           T      = Temperature moyenne journaliere / daily mean temperature (°C)
%           Lat    = Latitude de la station / station latitude (°)
%
%   Outputs
%       PetData.
%           T      = Mean air temperature
%           DL     = Maximum day light (h)
%
% Coded by A. Thiboult (2017)

%% Computation day of year (1-365)
doy = dayOfYear(Data.Date);

%% Computation
Lrad = (pi()*Data.Lat)/180; % (rad)
ds = 0.409*sin((2*pi()/365)*doy-1.39); % (rad)
omega = acos(-tan(Lrad).*tan(ds)); % (rad)
DL = 24/pi()*omega; % en h

%% Outputs
PetData.T      = Data.T;
PetData.DL     = DL;

end

