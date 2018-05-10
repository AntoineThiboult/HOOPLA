function [ETP] = Penman(Data)
% 
% [E] = Oudin(data)
% 
% Warning ! This function was not tested. It is likely that it is bugged 
%
% Computation of the potentatial evapotranspiration according to the Penman
% formula (1948)
% 
%   Inputs
%       Data.
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
%   Outputs
%       E = potential evapotranspiration (mm/3h)
% 
% Coded by G. Seiller
% Modified by A. Thiboult (2017)

Rn = Data.Rs.*(1-Data.alpha)-Data.sigma.*((Data.T+273.2).^4).*(0.56-0.25.*sqrt(Data.ed)).*(0.1+0.9*(Data.D./Data.DL))./8; % (MJ/m2/3h)
W = 2.62.*(1+0.537.*Data.U);
ETP = ((Data.delta./(Data.delta+Data.gamma)).*Rn+W.*(Data.gamma./(Data.delta+Data.gamma)).*(Data.es-Data.ed))./Data.lambda;
ETP = max(0,ETP);
end

