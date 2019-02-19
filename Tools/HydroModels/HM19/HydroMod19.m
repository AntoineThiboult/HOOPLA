function [Qsim, param, varargout] = HydroMod19( P, E, param )
%
% [Qsim, param, varargout] = HydroMod19( P, E, param )
%
% WAGENINGEN hydrological model (modified version)
%
% INPUTS (time series of daily observations [n,1])
% P 	= mean areal rainfall (mm)
% E 	= mean areal evapotranspiration (mm)
% Q 	= stream flow (mm)
% x 	= the ten model parameters (see "param" below) - [8,1]
%
% OUTPUTS
% Qs    = simulated stream flow (mm)
% perf  = model performances
% inter = WAGENINGEN's internal values
% param -->
% 	.x(1) = seuil de vidange des percolations
% 	.x(2) = capacité maximale du réservoir sol
% 	.x(3) = constante de vidange des infiltrations
% 	.x(4) = paramètre des remontées capillaires
% 	.x(5) = paramètre de dissociation des écoulements
% 	.x(6) = constante de vidange de l’écoulement rapide
% 	.x(7) = constante de vidange de l’écoulement lent
% 	.x(8) = délai
%
% 	.S = Réservoir sol (mm)
% 	.R = Réservoir de routage rapide (mm)
% 	.T = Réservoir de routage lent (mm)
%
% FOLLOWING
%  - Warmerdam, P.M., Kole, J., Chormanski, J., 1997. Modelling rainfall-runoff 
%    processes in the Hupselse beek research basin. In: IHP-V, Technical 
%    Documents in Hydrology. pp. 155–160. (Wageningen Agricultural University,
%    Netherlands)
%  - Perrin, C. (2000). Vers une amélioration d'un modèle global pluie-débit, 
%    PhD Thesis, Appendix 1, p. 459-462. Retrieved from 
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
DL = param.DL ;
HY = param.HY ;

%%% PRODUCTION PART

% Soil storage (S)
%
S = S+P;

if S >= x(1)
    Is = (S/x(2))*((S-x(1))/x(3));
    It = 0;
else
    Is = 0;
    It = (T/x(4))*(x(1)-S);
end

S = S+It-Is;

if S >= x(1)
    Es = E;
else
    Es = E*cos((pi/2)*((x(1)-S)/x(1)));
end

S = max(0,S-Es);

DIV = min(1,T/x(5));

T = T+(1-DIV)*Is;

R = R+DIV*Is;

%%% ROUTING
%
% Quick routing storage (R)
%
Qr = R/x(6);
R = R-Qr;

% Slow routing storage (T)
%
Qt = T/(x(6)*x(7));
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
    inter = [ S R T Is It Qr Qt ];
    
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