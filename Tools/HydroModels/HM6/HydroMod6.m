function [Qsim, param, varargout] = HydroMod6( P, E, param )
%
% [Qsim, param, varargout] = HydroMod6( P, E, param )
%
% HBV hydrological model (modified version)
%
% INPUTS (time series of daily observations [n,1])
% P 	= mean areal rainfall (mm)
% E 	= mean areal evapotranspiration (mm)
% Q 	= stream flow (mm)
% x 	= the nine model parameters (see "param" below) - [9,1]
%
% OUTPUTS
% Qs    = simulated stream flow (mm)
% perf  = model performances
% inter = HBV's internal values
% param -->
% 	.x(1) = Capacité du réservoir sol (mm)
% 	.x(2) = Seuil pour l'ETP
% 	.x(3) = Constante de vidange supérieure du réservoir intermédiaire
% 	.x(4) = Constante de vidange du réservoir souterrain
% 	.x(5) = Coefficient de percolation
% 	.x(6) = Constante de temps de l'hydrogramme triangulaire
% 	.x(7) = Exposant B
% 	.x(8) = Seuil d'écoulement du réservoir intermédiaire
% 	.x(9) = Constante de vidange inférieure du réservoir intermédiaire
%
% 	.S = Réservoir de sol (mm)
% 	.R = Réservoir intermédiaire (mm)
% 	.T = Réservoir souterrain (mm)
%
% 	.UH  = Unit hydrograph
% 	.H   = Hydrograph values (mm) - updated at each time step
%
% FOLLOWING
%  - Bergström, S., Forsman, A., 1973. Development of a conceptual 
%    deterministic rainfall-runoff model. Nord. Hydrol. 4, 147–170. 
%    (Swedish Meteorological and Hydrological Institute, Norrköping,
%    Sweden)
%  - Perrin, C. (2000). Vers une amélioration d'un modèle global pluie-débit, 
%    PhD Thesis, Appendix 1, p. 366-371. Retrieved from 
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
UH = param.UH ;
H  = param.H ;

%%% PRODUCTION PART

% Soil storage (S) routine (5 steps)
%
Pr = 0;

P5 = P/5;
E5 = E/5;

Pr1 = P5*((min(1,S/x(1)))^x(7));
Pr = Pr + Pr1;
S = S + (P5 - Pr1);
Es1 = min(S,E5*(S/x(2)));
S = S-Es1;

Pr2 = P5*((min(1,S/x(1)))^x(7));
Pr = Pr + Pr2;
S = S + (P5 - Pr2);
Es2 = min(S,E5*(S/x(2)));
S = S-Es2;

Pr3 = P5*((min(1,S/x(1)))^x(7));
Pr = Pr + Pr3;
S = S + (P5 - Pr3);
Es3 = min(S,E5*(S/x(2)));
S = S-Es3;

Pr4 = P5*((min(1,S/x(1)))^x(7));
Pr = Pr + Pr4;
S = S + (P5 - Pr4);
Es4 = min(S,E5*(S/x(2)));
S = S-Es4;

Pr5 = P5*((min(1,S/x(1)))^x(7));
Pr = Pr + Pr5;
S = S + (P5 - Pr5);
Es5 = min(S,E5*(S/x(2)));
S = S-Es5;

%%% ROUTING
%
% Routing storage (R)
%
R = R+Pr;
Qr1 = max(0,(R-x(8))/x(3));
R = R-Qr1;
Qr2 = R/(x(3)*x(9));
R = R-Qr2;
Ir = min(R,x(5));
R = R-Ir;

% Groundwater storage (T)
%
T = T+Ir;
Qt = T/x(4);
T = T-Qt;

% Total discharge
%
Q = Qr1 + Qr2 + Qt;

% Mise à jour de l'hydrogramme
%
H = [H(2:end); 0] + UH * Q ;

Qsim = max( [0; H(1)] ) ;

param.S = S;
param.R = R;
param.T = T;
param.H = H ;

if nargout>2
    inter = [ S R T Pr Qr1 Qr2 Qt ];
    
    % Flow components
    qsf = 0 ;
    qs1 = 0 ; qs2 = 0 ; qs = qs1+qs2 ;
    qrs1 = 0 ; qrs2 = 0 ; qrs = qrs1+qrs2 ;
    qss1 = Qr1+Qr2 ; qss2 = 0 ; qss = qss1+qss2 ;
    qrss1 = 0 ; qrss2 = 0 ; qrss = qrss1+qrss2 ;
    qn = Qt ;
    qr = 0 ;
    
    interq = [ qsf qs1 qs2 qs qrs1 qrs2 qrs qss1 qss2 qss qrss1 qrss2 qrss qn qr ];
    
    % Outputs
    varargout{1}=inter;
    varargout{2}=interq;
end
