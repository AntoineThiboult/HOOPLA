function [Qsim, param, varargout] = HydroMod9( P, E, param )
%
% [Qsim, param, varargout] = HydroMod9( P, E, param )
%
% MARTINE hydrological model (modified version)
%
% INPUTS (time series of daily observations [n,1])
% P       = mean areal rainfall (mm)
% E       = mean areal evapotranspiration (mm)
% Q       = stream flow (mm)
% x       = the seven model parameters (see "param" below) - [7,1]
%
% OUTPUTS
% Qs      = simulated stream flow (mm)
% perf    = model performances
% inter   = MARTINE's internal values
% param -->
%   .x(1) = capacité du réservoir superficiel
%   .x(2) = capacité du réservoir intermédiaire
%   .x(3) = capacité du réservoir de routage quadratique
%   .x(4) = constante de vidange du réservoir souterrain
%   .x(5) = coefficient de partage
%   .x(6) = délai
%   .x(7) = constante de vidange du réservoir intermédiaire
%   .S    = Réservoir de superficiel
%   .T    = Réservoir de intermédiaire
%   .L    = Réservoir souterrain
%   .R    = Réservoir de routage quadratique
%
% FOLLOWING
% Mazenc et al. (1984)
% Bureau de Recherches Géologiques et Minières, Orléans, France
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

%% Modeling
% Get parameters
%
x   = param.x ;
S   = param.S ;
T   = param.T ;
L   = param.L ;
R   = param.R ;
DL = param.DL ;
HY = param.HY ;

%%%
%%% PRODUCTION AND ROUTING
%%%

% Soil storage (S)
%
S = S+P;
Pr = max(0,S-x(1));
S = S-Pr;
Es = min(S,E);
S = S-Es;
Et = E-Es;

% Direct routing storage (R)
%
R = R+x(5)*Pr;
Qr = (R^2)/(R+x(3));
R = R-Qr;

% Lower storage (T)
%
T = max(0,T-Et);
T = T+(1-x(5))*Pr;
Qt1 = max(0,T/x(7));
T = T-Qt1;
Qt2 = max(0,T-x(2));
T = T-Qt2;

% Groundwater storage (L)
%
L = L+Qt1+Qt2;
Ql = L/x(4);
L = L-Ql;

% Total discharge
%
HY = [HY(2:end); 0] + DL * (Ql+Qr) ;
Qsim = max( [0; HY(1)] ) ;

% Data
param.S  = S ;
param.T  = T ;
param.L  = L ;
param.R  = R ;
param.HY = HY ;

if nargout>2
    inter = [ S T L R Pr Es Et Qt1 Qt2 Ql Qr ] ;
    
    % Flow components
    qsf = 0 ;
    qs1 = 0 ; qs2 = 0 ; qs = qs1+qs2 ;
    qrs1 = 0 ; qrs2 = 0 ; qrs = qrs1+qrs2 ;
    qss1 = 0 ; qss2 = 0 ; qss = qss1+qss2 ;
    qrss1 = 0 ; qrss2 = 0 ; qrss = qrss1+qrss2 ;
    qn = Ql ;
    qr = Qr ;
    
    interq = [ qsf qs1 qs2 qs qrs1 qrs2 qrs qss1 qss2 qss qrss1 qrss2 qrss qn qr ];
    
    % Outputs
    varargout{1}=inter;
    varargout{2}=interq;
end