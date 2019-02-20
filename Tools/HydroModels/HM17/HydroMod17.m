function [Qsim, param, varargout] = HydroMod17( P, E, param )
% 
% [Qsim, param, varargout] = HydroMod17( P, E, param )
%
% TANK hydrological model (modified version)
%
% INPUTS
% P 	= mean areal rainfall (mm)
% E 	= mean areal evapotranspiration (mm)
% param	= the seven model parameters (see "param" below) - [7,1]
%
% OUTPUTS
% Qsim    = simulated stream flow (mm)
% inter   = HydroMod17's internal values (varargout 1)
% interq  = HydroMod17's internal flow components (varargout 2)
% param -->
% 	.x(1) = Upper flow threshold of soil reservoir
% 	.x(2) = Lower flow threshold of soil reservoir
% 	.x(3) = Fast emptying constant
% 	.x(4) = Intermediate emptying constant
% 	.x(5) = Delay
% 	.x(6) = PET correction coefficient
% 	.x(7) = Slow emptying constant
%
% 	.S = Surface reservoir state (mm)
% 	.R = Upper soil reservoir state (mm)
% 	.T = Lower soil reservoir state (mm)
% 	.L = Ground reservoir state (mm)
%
% FOLLOWING
%  - Sugawara, M., 1979. Automatic calibration of the tank model. 
%    Hydrol. Sci. 24, 375–388. (National Research Centre for Disaster 
%    Prevention, Tokyo, Japan)
%  - Perrin, C. (2000). Vers une amélioration d'un modèle global pluie-débit, 
%    PhD Thesis, Appendix 1, p. 443-447. Retrieved from 
%    https://tel.archives-ouvertes.fr/tel-00006216
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
L = param.L ;
DL = param.DL ;
HY = param.HY ;

% Surface storage (S)
%
S = S+P;
E1 = E*x(6);
Qs1 = max(0,(S-(x(1)+x(2)))/x(3));
S = S-Qs1;
Qs2 = max(0,(S-x(2))/x(3));
S = S-Qs2;
Is = S/x(3);
S = S-Is;
Es = min(E1,S);
S = S-Es;
E2 = E1-Es;

% Upper soil storage (R)
%
R = R+Is;
Qr = max(0,(R-x(2))/(x(3)*x(4)));
R = R-Qr;
Ir = R/(x(3)*x(4));
R = R-Ir;
Er = min(E2,R);
R = R-Er;
E3 = E2-Er;

% Lower soil storage (T)
%
T = T+Ir;
Qt = max(0,(T-x(2))/(x(3)*x(4)*x(7)));
T = T-Qt;
It = T/(x(3)*x(4)*x(7));
T = T-It;
Et = min(E3,T);
T = T-Et;
E4 = E3-Et;

% Groundwater storage (L)
%
L = L+It;
Ql = L/(x(3)*x(4)*(x(7)^2)); % x(7)^2 sur graph et que x(7) sur formules
L = L-Ql;
El = min(E4,L);
L = L-El;

% Total discharge
%
HY = [HY(2:end); 0] + DL * (Qs1 + Qs2 + Qr + Qt + Ql) ;
Qsim = max( [0; HY(1)] ) ;

param.S = S;
param.R = R;
param.T = T;
param.L = L;
param.HY = HY ;

if nargout>2
    inter = [ S R T L Is Ir It Qs1 Qs2 Qr Qt Ql ];
    
    % Flow components
    qsf = 0 ;
    qs1 = Qs1+Qs2 ; qs2 = 0 ; qs = qs1+qs2 ;
    qrs1 = 0 ; qrs2 = 0 ; qrs = qrs1+qrs2 ;
    qss1 = Qr ; qss2 = Qt ; qss = qss1+qss2 ;
    qrss1 = 0 ; qrss2 = 0 ; qrss = qrss1+qrss2 ;
    qn = Ql ;
    qr = 0 ;
    
    interq = [ qsf qs1 qs2 qs qrs1 qrs2 qrs qss1 qss2 qss qrss1 qrss2 qrss qn qr ];
    
    % Outputs
    varargout{1}=inter;
    varargout{2}=interq;
end