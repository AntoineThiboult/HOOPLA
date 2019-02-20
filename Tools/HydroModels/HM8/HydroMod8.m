function [Qsim, param, varargout] = HydroMod8( P, E, param )
%
% [Qsim, param, varargout] = HydroMod8( P, E, param )
%
% IHACRES hydrological model (modified version)
%
% INPUTS
% P 	= mean areal rainfall (mm)
% E 	= mean areal evapotranspiration (mm)
% param = the seven model parameters (see "param" below) - [7,1]
%
% OUTPUTS
% Qsim    = simulated stream flow (mm)
% inter   = HydroMod8's internal values (varargout 1)
% interq  = HydroMod8's internal flow components (varargout 2)
% param -->
% 	.x(1) = Forcing parameter 1/C (computed to ensure water balance over the whole period)
% 	.x(2) = Flow distribution paramter
% 	.x(3) = Emptying constant of the fast routing reservoir
% 	.x(4) = Emptying constant of the slow routing reservoir (X3.X4)
% 	.x(5) = Delay
% 	.x(6) = f parameter (temperature modulation factor)
% 	.x(7) = Tw parameter (characteristic watershed drying constant)
%
% 	.S = Soil moisture reservoir state (mm)
% 	.T = Slow routing reservoir state (mm)
% 	.R = Fast routing reservoir state (mm)
%
% FOLLOWING
%  - Littlewood et al. (1997), The PC version of IHACRES for catchment-scale 
%    rainfall-streamflow modelling. Version 1.0. User Guide. Institute of 
%    Hydrology(Ed.), 89 p
%  - Perrin, C. (2000). Vers une amélioration d'un modèle global pluie-débit, 
%    PhD Thesis, Appendix 1, p. 376-382. Retrieved from 
%    https://tel.archives-ouvertes.fr/tel-00006216
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

%% Modeling
% Get parameters
%
x  = param.x ;
S  = param.S ;
T  = param.T ;
R  = param.R ;
DL = param.DL ;
HY = param.HY ;

%%%
%%% PRODUCTION PART
%%%

% Soil moisture accounting
%
XS = S;
El = max(0,x(7)-(E/x(6)));
S = XS+(P/x(1))-(XS/exp(El));

% Percolation
%
Pr = 0.5*(XS+S)*P;

%%%
%%% ROUTING
%%%

% Quick routing

R = R+x(2)*Pr;
Qr = R/x(3);
R = R-Qr;

% Slow routing

T = T+(1-x(2))*Pr;
Qt = T/(x(3)*x(4));
T = T-Qt;

% Total discharge

HY = [HY(2:end); 0] + DL * (Qt + Qr) ;
Qsim = max( [0; HY(1)] ) ;


param.S  = S ;
param.R  = R ;
param.T  = T ;
param.HY = HY ;

if nargout>2
    inter = [ XS El S Pr T R Qt Qr ] ;
    
    % Flow components
    qsf = 0 ;
    qs1 = 0 ; qs2 = 0 ; qs = qs1+qs2 ;
    qrs1 = 0 ; qrs2 = 0 ; qrs = qrs1+qrs2 ;
    qss1 = 0 ; qss2 = 0 ; qss = qss1+qss2 ;
    qrss1 = 0 ; qrss2 = 0 ; qrss = qrss1+qrss2 ;
    qn = Qt ;
    qr = Qr ;
    
    interq = [ qsf qs1 qs2 qs qrs1 qrs2 qrs qss1 qss2 qss qrss1 qrss2 qrss qn qr ];
    
    % Outputs
    varargout{1}=inter;
    varargout{2}=interq;
end