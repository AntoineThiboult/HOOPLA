function [PetData] = ini_Penman(~, Data)
% [PetData] = ini_Penman(Data)
%
% Warning ! This function was not tested. It is likely that it is bugged 
%
% Compute the prerequisite values that are necessary for the computation of Oudin
%
%   Inputs:
%       Data.
%           Date   = Date / date (datenum)
%           Glorad = Rayonnement solaire incident / measured incoming solar radiation (W/m2)
%           Lat    = Latitude de la station / station latitude (°)
%           Relhum = Humidité relative / relative humidity (%)
%           Tmax   = Température max journalière / daily max temperature (°C)
%           T       = Température moyenne journalière / daily mean temperature (°C)
%           Tmin   = Température min journalière / daily min temperature (°C)
%           Wndspd = Vitesse du vent à 2 mètres / wind speed at 2 meters (m/s)
%           z      = Altitude de la station / station altitude (m)
%   Outputs:
%       PetData.
%           Rs      = Measured incoming solar radiation (MJ/m2/j)
%           alpha   = Albedo
%           sigma   = Stefan-Boltzmann constant
%           T       = Mean air temperature
%           ed      = Dew point temperature
%           D       = Effective day light
%           DL      = Maximum day light (h)
%           U       = Wind speed at 2 meters (m/s)
%           delta   = slope of the saturating vapor pressure curve (kPa)
%           gamma   = Psychometric constant (kPa/degCel)
%           es      = e(Tmax)+e(Tmin/2) Saturating vapor pressure (kPa)
%           lambda = Latent vaporization energy (MJ/kg)
%
% Coded by A. Thiboult (2017)

%% Computation day of year (1-365)
doy = dayOfYear(Data.Date);

%% Constants
Gsc = 0.082; % in MJ/m2/min
lambda = 2.501-0.002361.*Data.T; % in MJ/kg
P = 101.3*(((293-0.0065*Data.z)/293)^5.26); % in kPa
gamma = 0.0016286*(P./lambda); % in kPa/°C
sigma = 4.903*10^(-9) ; % in MJ/m2/j
rs = 70; % in s/m (70 for crops of grass 0.12m and 45 for alfalfa 0.5 m)

%% Calculs préalables
Rs = Data.Glorad * 0.0864; % de W/m2 à MJ/m2/day
alpha = 0.29+0.06.*sin((doy+96)./57.3);
Td = (237.3*(((17.27.*Data.T)./(237.3+Data.T))+log(Data.Relhum/100)))./(17.27-(((17.27.*Data.T)./(237.3+Data.T))+log(Data.Relhum/100))); % in °C Formula de Tetens
ea = 0.6108*exp((17.27.*Data.T)./(237.3+Data.T)) ; % in kPa
ed = 0.6108*exp((17.27.*Td)./(237.3+Td)) ; % in kPa
es = ((0.6108*exp((17.27.*Data.Tmax)./(237.3+Data.Tmax)))+(0.6108*exp((17.27.*Data.Tmin)./(237.3+Data.Tmin))))/2; % Formula ASCE
delta = (4098.*ea)./((237.3+Data.T).^2) ; % in kPa
Lrad = (pi()*Data.Lat)/180; % in rad
ds = 0.409*sin((2*pi()/365)*doy-1.39); % in rad
dr = 1+0.033*cos(doy*2*pi()/365); % without unit
omega = acos(-tan(Lrad).*tan(ds)); % in rad
Re = (24*60/pi())*Gsc.*dr.*(omega*sin(Lrad).*sin(ds)+cos(Lrad).*cos(ds).*sin(omega)); % in MJ/m2/j
DL = 24/pi()*omega; % in h
D = ((Rs./Re)-0.25).*(1/0.50).*DL; % Angstrom formula with as = 0.25 and bs = 0.50

%% Outputs
PetData.rs = rs;
PetData.alpha = alpha;
PetData.sigma = sigma;
PetData.T = Data.T;
PetData.ed = ed;
PetData.D = D;
PetData.DL = DL;
PetData.U = Data.Wndspd;
PetData.delta = delta;
PetData.gamma = gamma;
PetData.es = es;
PetData.lambda = lambda;
end