function [Qsim, param, varargout] = HydroMod5( P, E, param )
%
% [Qsim, param, varargout] = HydroMod5( P, E, param )
% 
% GR4H hydrological model
%
% INPUTS (time series of daily observations [n,1])
% P       = mean areal rainfall (mm)
% E       = mean areal evapotranspiration (mm)
% Q       = stream flow (mm)
% x       = the four model parameters (see "param" below) - [4,1]
%
% OUTPUTS
% Qsim    = simulated stream flow (mm)
% inter   = GR4J's internal values
% interq  = GR4J's flow components
% param -->
%   .x(1) = maximum capacity of the production store (mm)
%   .x(2) = groundwater exchange coefficient (mm)
%   .x(3) = one-day-ahead maximum capacity of the routing store (mm)
%   .x(4) = time base of unit hydrograph UH1 (day)
%   .B    = fraction of Pr routed by UH1 (fixed to 0.9)
%   .S    = production store level (mm)
%   .R    = routing store level (mm)
%   .UH1  = Unit hydrograph 1 - rapid flow
%   .UH2  = Unit hydrograph 2 - slower flow
%   .H1   = Hydrograph 1 values (mm) - updated at each time step
%   .H2   = Hydrograph 2 values (mm) - updated at each time step
%
% FOLLOWING
% Perrin C, Michel C & Andréassian V. 2003. Improvement of a parsimonious
%    model for streamflow simulation. J. Hydrol. 279, 275-289.
%
% Programmed by F. Anctil, Univ. Laval (01-2009)
% Modified by G. Seiller (2013),
% Modified by A. Thiboult (2016)

%% Modeling
% get parameters
%
x   = param.x ;
S   = param.S ;
R   = param.R ;
B   = param.B ;
UH1 = param.UH1 ;
UH2 = param.UH2 ;
H1  = param.H1 ;
H2  = param.H2 ;

%%%
%%% PRODUCTION PART
%%%

% Net inputs
%
if P >= E
    Pn = P - E ;
    En = 0 ;
else
    En = E - P ;
    Pn = 0 ;
end

% Soil moisture accounting
%
if Pn > 0
    tilap1 = S / x(1) ;
    tilap2 = tanh( Pn / x(1) ) ;
    Ps = ( x(1)*( 1-tilap1.^2 )*tilap2 ) / ( 1+tilap1*tilap2 ) ;
    Es = 0 ;
else
    tilap1 = S / x(1) ;
    tilap2 = tanh( En / x(1) ) ;
    Es = ( S*(2-tilap1)*tilap2 ) / ( 1+(1-tilap1)*tilap2 ) ;
    Ps = 0 ;
end
S = S - Es + Ps ;

% Percolation
%
if x(1)/S > 0.001
    Perc =  S * ( 1- ( 1+ (4*S/21/x(1)).^4).^-0.25 ) ;
    S = S - Perc ;
else
    Perc = 0 ;
end
Pr = Pn - Ps + Perc ;

%%%
%%% ROUTING
%%%

% Mise à jour des hydrogrammes 1 & 2
%
H1 = [H1(2:end); 0] + UH1 * Pr * B ;
H2 = [H2(2:end); 0] + UH2 * Pr * (1-B) ;

% Calcul de l'échange
%
F = x(2) * (R/x(3))^3.5 ;

% Mise à jour du niveau du réservoir de routage
%
R = max( [0.001*x(3); R + H1(1) + F] ) ;

% Calcul de la vidange du réservoir de routage et mise à jour du niveau
%
Qr = R * ( 1- ( 1+ (R/x(3)).^4 ).^-0.25 ) ;
R = R - Qr ;

% Calcul de la composante pseudo-direct de l'écoulement
%
Qd = max( [0; H2(1) + F] ) ;

% Calcul du débit total
%
Qsim = Qr + Qd ;

param.S  = S ;
param.R  = R ;
param.H1 = H1 ;
param.H2 = H2 ;

if nargout>2
    inter = [ Es Ps S Perc Pr F R Qr Qd ] ;
    
    % Flow components
    qsf = 0 ;
    qs1 = Qd ; qs2 = 0 ; qs = qs1+qs2 ;
    qrs1 = 0 ; qrs2 = 0 ; qrs = qrs1+qrs2 ;
    qss1 = 0 ; qss2 = 0 ; qss = qss1+qss2 ;
    qrss1 = 0 ; qrss2 = 0 ; qrss = qrss1+qrss2 ;
    qn = 0 ;
    qr = Qr ;
    
    interq = [ qsf qs1 qs2 qs qrs1 qrs2 qrs qss1 qss2 qss qrss1 qrss2 qrss qn qr ];
    
    % Outputs 
    varargout{1}=inter;
    varargout{2}=interq;    
end