function [Qsim, param, varargout] = HydroMod4( P, E, param )
%
% [Qsim, param, varargout] = HydroMod4( P, E, param )
%
% GARDENIA hydrological model (modified version)
%
% INPUTS (time series of daily observations [n,1])
% P 	= mean areal rainfall (mm)
% E 	= mean areal evapotranspiration (mm)
% Q 	= stream flow (mm)
% x 	= the six model parameters (see "param" below) - [6,1]
%
% OUTPUTS
% Qs    = simulated stream flow (mm)
% perf  = model performances
% inter = GARDENIA's internal values
% param -->
% 	.x(1) = capacité du réservoir de surface
% 	.x(2) = constante de percolations linéaires
% 	.x(3) = paramètre de vidange latérale du réservoir sol
% 	.x(4) = constante de vidange linéaire du réservoir souterrain
% 	.x(5) = coefficient de correction des ETP
% 	.x(6) = délai
%
% 	.S = Réservoir de surface (mm)
% 	.R = Réservoir de sol (mm)
% 	.T = Réservoir souterrain (mm)
%
% FOLLOWING
% Thiery (1982)
% Bureau de Recherches Géologiques et Minières (BRGM, ), Orléans (FR)
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

%% Modeling
% Get parameters
%
x = param.x ;
S = param.S ;
R = param.R ;
T = param.T ;
DL = param.DL ;
HY = param.HY ;

%%%
%%% PRODUCTION PART
%%%

% Surface storage (S)
%
S = S+P;
Pr = max(0,S-x(1));
S = S-Pr;
Es = x(5)*E;
S = max(0,S-Es);

%%%
%%% ROUTING
%%%

% Soil storage (R)
%
R = R+Pr;
Qr = (R^2)/(R+x(2)*x(3));
R = R-Qr;
Ir = R/x(2);
R = R-Ir;

% Groundwater storage (T)
%
T = T+Ir;
Qt = T/x(4);
T = T-Qt;

% Total discharge
%
HY = [HY(2:end); 0] + DL * (Qt+Qr) ;
Qsim = max( [0; HY(1)] ) ;

param.S = S;
param.R = R;
param.T = T;
param.HY = HY ;

if nargout>2
    inter = [ S R T Es Pr Ir Qr Qt ];
    
    % Flow components
    qsf = 0 ;
    qs1 = Qr ; qs2 = 0 ; qs = qs1+qs2 ;
    qrs1 = 0 ; qrs2 = 0 ; qrs = qrs1+qrs2 ;
    qss1 = 0 ; qss2 = 0 ; qss = qss1+qss2 ;
    qrss1 = 0 ; qrss2 = 0 ; qrss = qrss1+qrss2 ;
    qn = Qt ;
    qr = 0 ;
    
    interq = [ qsf qs1 qs2 qs qrs1 qrs2 qrs qss1 qss2 qss qrss1 qrss2 qrss qn qr ];
    
    % Outputs
    varargout{1}=inter;
    varargout{2}=interq;
end