function [Qsim, param, varargout] = HydroMod2( P, E, param )
%
% [Qsim, param, varargout] = HydroMod2( P, E, param )
%
% CEQUEAU hydrological model (modified version)
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
% inter = CEQUEAU's internal values
% param -->
% 	.x(1) = seuil d'infiltration (mm)
% 	.x(2) = seuil de vidange du premier réservoir (mm)
% 	.x(3) = constante de vidange d’infiltration
% 	.x(4) = constante de vidange latérale supérieure du réservoir sol
% 	.x(5) = capacité maximum du réservoir souterrain (mm)
% 	.x(6) = délai
% 	.x(7) = seuil de vidange du réservoir souterrain (mm)
% 	.x(8) = constante de vidange latérale inférieure du réservoir sol
% 	.x(9) = constante de vidange inférieure du réservoir souterrain
%
% 	.S = Réservoir de surface (mm)
% 	.T = Réservoir souterrain (mm)
%
% FOLLOWING
% INRS Eau
% Girard et al. (1972)
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

%% Modeling
% get parameters
%
x = param.x ;
S = param.S ;
T = param.T ;
DL = param.DL ;
HY = param.HY ;

%%% PRODUCTION PART

% Net inputs
%
S = S+P;
Es = min(S,E*min(1,(2*S)/x(5)));
S = S-Es;
Er = E-Es;

% Percolation
%
Is = max(0,S-x(1))/x(3);
S = S-Is;

% Surface routing and updates
%
Qs2 = max(0,S-x(2))/x(4);
S = S-Qs2;
Qs3 = S/(x(4)*x(8));
S = S-Qs3;
Qs1 = max(0,S-x(5));
S = S-Qs1;

% Groundwater routing and updates
%
T = T+Is;
Qt1 = max(0,T-x(7))/(x(4)*x(9));
T = T-Qt1;
Qt2 = T/(x(4)*x(8)*(x(9)^2));
T = T-Qt2;
Et = min(T,Er*min(1,T/x(7)));
T = T-Et;

% Total discharge
%
HY = [HY(2:end); 0] + DL * (Qs1+Qs2+Qs3+Qt1+Qt2) ;
Qsim = max( [0; HY(1)] ) ;

param.S = S;
param.T = T;
param.HY = HY ;

if nargout>2
    inter = [ P E S T Es Er Et Is Qs1 Qs2 Qs3 Qt1 Qt2 ];
    
    % Flow components
    qsf = 0 ;
    qs1 = Qs1+Qs2+Qs3 ; qs2 = 0 ; qs = qs1+qs2 ;
    qrs1 = 0 ; qrs2 = 0 ; qrs = qrs1+qrs2 ;
    qss1 = 0 ; qss2 = 0 ; qss = qss1+qss2 ;
    qrss1 = 0 ; qrss2 = 0 ; qrss = qrss1+qrss2 ;
    qn = Qt1+Qt2 ;
    qr = 0 ;
    
    interq = [ qsf qs1 qs2 qs qrs1 qrs2 qrs qss1 qss2 qss qrss1 qrss2 qrss qn qr ];
    
    % Outputs
    varargout{1}=inter;
    varargout{2}=interq;
end
