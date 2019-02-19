function [Qsim, param, varargout] = HydroMod3( P, E, param )
%
% [Qsim, param, varargout] = HydroMod3( P, E, param )
%
% CREC hydrological model (modified version)
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
% inter = CREC's internal values
% param -->
% 	.x(1) = constante de vidange du réservoir souterrain
% 	.x(2) = paramètre de percolation linéaire du réservoir sol
% 	.x(3) = paramètre de séparation de la pluie brute
% 	.x(4) = paramètre de séparation de la pluie brute et de rendement d'ETP
% 	.x(5) = paramètre de vidange linéaire du réservoir sol
% 	.x(6) = délai
%
% 	.S = Réservoir de sol (mm)
% 	.R = Réservoir de routage de sol (mm)
% 	.T = Réservoir souterrain (mm)
%
% FOLLOWING
%  - Cormary, Y., Guilbot, A., 1973. Étude des relations pluie-débit sur 
%    trois bassins versants d’investigation. In: publications, I. (Ed.), 
%    Design of water resources projects with inadequate data: IAHS Madrid 
%    Symposium. pp. 265–279. (Laboratoire d’Hydrologie Mathématique, 
%    Université des Sciences, Montpellier, France.)
%  - Perrin, C. (2000). Vers une amélioration d'un modèle global pluie-débit, 
%    PhD Thesis, Appendix 1, p. 327-332. Retrieved from 
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
XF = param.XF;
DL = param.DL ;
HY = param.HY ;

%%%
%%% PRODUCTION PART
%%%

% Net inputs
%
Pr = P/(1+exp((x(3)-S)/x(4)));
Ps = P-Pr;

% Soil storage (S)
%
S = S+Ps;
Es = E*(1-exp(-(S/XF)));
S = max(0,S-Es);

%%%
%%% ROUTING
%%%

% Surface storage (R)
%
R = R+Pr;
Qr = (R^2)/(R+x(1));
R = R-Qr;
Ir = R/x(5);
R = R-Ir;

% Groundwater storage (T)
%
T = T+Ir;
Qt = T/x(2);
T = T-Qt;

% Total discharge
%
HY = [HY(2:end); 0] + DL * (Qr+Qt) ;
Qsim = max( [0; HY(1)] ) ;

param.S = S;
param.R = R;
param.T = T;
param.HY = HY ;

if nargout>2
    inter = [ S R T Ps Pr Es Ir Qr Qt];
    
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