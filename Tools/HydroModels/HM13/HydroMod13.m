function [Qsim, param, varargout] = HydroMod13( P, E, param )
%
% [Qsim, param, varargout] = HydroMod13( P, E, param )
%
% PDM8 hydrological model (modified version)
%
% INPUTS
% P 	= mean areal rainfall (mm)
% E 	= mean areal evapotranspiration (mm)
% param	= the eight model parameters (see "param" below) - [8,1]
%
% OUTPUTS
% Qsim    = simulated stream flow (mm)
% inter   = HydroMod13's internal values (varargout 1)
% interq  = HydroMod13's internal flow components (varargout 2)
% param -->
% 	.x(1) = Maximum capacity of production reservoir
% 	.x(2) = Spatial variability coefficient of storage capacity (Bexp)
% 	.x(3) = Drainage parameter (Alpha)
% 	.x(4) = Delay
% 	.x(5) = Emptying constant of the ground reservoir
% 	.x(6) = Emptying constant of series routing reservoirs (Rq)
% 	.x(7) = Rainfall correction factor
% 	.x(8) = Drainage constant
%
% 	.S = Soil reservoir state (mm)
% 	.T = Ground reservoir state (mm)
% 	.M = Routing reservoir in series 1 (mm)
% 	.N = Routing reservoir in series 2 (mm)
%
% FOLLOWING
%  - Moore, R.J., Clarke, R.T., 1981. A distribution function approach to 
%    rainfall runoff modeling. Water Resour. Res. 17, 1367–1382.
%    (Institute of Hydrology, Wallingford, UK)
%  - Perrin, C. (2000). Vers une amélioration d'un modèle global pluie-débit, 
%    PhD Thesis, Appendix 1, p. 420-424. Retrieved from 
%    https://tel.archives-ouvertes.fr/tel-00006216
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

%% Modeling
% Get parameters
%
x = param.x ;
S = param.S ;
T = param.T ;
M = param.M ;
N = param.N ;
DL = param.DL ;
HY = param.HY ;

%%%
%%% PRODUCTION PART
%%%


Cmax = x(1);
Bexp = x(2);
Alpha = x(3);
Rq = x(6);

% Net inputs
%

% Correction de la pluie
P1 = P*x(7);
Xs = S;

Ctprev = Cmax*(1-(1-(Bexp+1)*S/Cmax)^(1/(Bexp+1)));
Ut1 = max(0,P1-Cmax+Ctprev);
Pn = P1-Ut1;
Dum = min(1,(Ctprev+Pn)/Cmax);
S = Cmax/(Bexp+1)*(1-(1-Dum)^(Bexp+1));
Ut2 = max(0,Pn-(S-Xs));

% Evaporation
S = max(0,S-E*(1-(1-S/Cmax*(Bexp+1))^2));

% Drainage
if S > Cmax/(Bexp+1)*Alpha
    Drg = (S-Cmax/(Bexp+1)*Alpha)/x(8);
else
    Drg = 0;
end

S = S-Drg;
Uq = Ut2+Ut1;
Us = Drg;

%%%
%%% ROUTING
%%%

% Routage rapide
M = M+Uq;
Q1 = M/Rq;
M = M-Q1;

N = N+Q1;
Q2 = N/Rq;
N = N-Q2;

% Routage lent (réservoir cubique)
T = T+Us;
Qt = T*(1-(1+(T/x(5))^2)^(-1/2));
T = T-Qt;

% Total discharge
%
HY = [HY(2:end); 0] + DL * (Qt+Q2) ;
Qsim = max( [0; HY(1)] ) ;

param.S = S;
param.T = T;
param.M = M;
param.N = N;
param.HY = HY ;

if nargout>2
    inter = [ S T M N Drg Q1 Q2 Qt ];
    
    % Flow components
    qsf = 0 ;
    qs1 = 0 ; qs2 = 0 ; qs = qs1+qs2 ;
    qrs1 = 0 ; qrs2 = Q2 ; qrs = qrs1+qrs2 ;
    qss1 = 0 ; qss2 = 0 ; qss = qss1+qss2 ;
    qrss1 = 0 ; qrss2 = 0 ; qrss = qrss1+qrss2 ;
    qn = Qt ;
    qr = 0 ;
    
    interq = [ qsf qs1 qs2 qs qrs1 qrs2 qrs qss1 qss2 qss qrss1 qrss2 qrss qn qr ];
    
    % Outputs
    varargout{1}=inter;
    varargout{2}=interq;
end