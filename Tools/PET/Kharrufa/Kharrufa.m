function [E] = Kharrufa(Data)
% 
% [E] = Kharrufa(data)
%  
%   Inputs
%       Data.
%           T      = Mean air temperature
%           DL     = Maximum day light (h)
%
%   Outputs
%       E = potential evapotranspiration (mm/3h)
%
% Computation of the potentatial evapotranspiration according to the
% Kharrufa formula 
% 
% Coded by G. Seiller
% Modified by A. Thiboult (2017)

p = 100*Data.DL/(365*12*8); % (3h)
k = 0.34; % Kharrufa empirical coeff 
E = k.*p.*(max(0,Data.T).^(1.3)) ;% (mm/3h)
E = max(0,E);
end