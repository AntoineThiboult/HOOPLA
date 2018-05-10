function [E] = Oudin(Data)
% 
% [E] = Oudin(data)
% 
% Computation of the potentatial evapotranspiration according to the Oudin
% formula 
% 
%   Inputs
%       Data.
%           Re     = Extraterrestrial radiation (MJ/m2/3h)
%           lambda = Latent vaporization energy (MJ/kg)
%           rho    = Water volumic mass (= 1000 kg/L)
%           T      = Mean air temperature
%           DL     = Maximum day light (h)
%
%   Outputs
%       E = potential evapotranspiration (mm/3h)
%
% Reference : Which potential evapotranspiration input for a lumped 
%   rainfall-runoff model? Part 2 - Towards a simple and efficient 
%   potential evapotranspiration model for rainfall-runoff modelling,
%   Oudin et al., Journal of Hydrology, 2005, 290-306, 303
% 
% Coded by G. Seiller
% Modified by A. Thiboult (2017)

E = (Data.Re./(Data.lambda.*Data.rho)).*((Data.T+5)/100) ; % (m/j)
E = E*1000 ; % (mm/j)
E = max(0,E);

end
