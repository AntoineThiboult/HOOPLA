function [Qsim, param, varargout] = HydroMod7( P, E, param )
%
% [Qsim, param, varargout] = HydroMod7( P, E, param )
%
% HYMOD hydrological model (modified version)
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
% inter = HYMOD's internal values
% param -->
% 	.x(1) = Capacité maximale du réservoir de sol (Cmax)
% 	.x(2) = Variabilité spatiale de la capacité d'humidité du sol (Bexp)
% 	.x(3) = Facteur de répartition des écoulements rapide/lent (alpha)
% 	.x(4) = délai
% 	.x(5) = Constante de vidange du réservoir de routage lent (Rs)
% 	.x(6) = Constante de vidange des réservoirs de routage rapide (Rq)
%
% 	.S = Réservoir de sol (mm)
% 	.R1 = Réservoir de routage rapide 1 (mm)
% 	.R2 = Réservoir de routage rapide 2 (mm)
% 	.R3 = Réservoir de routage rapide 2 (mm)
% 	.T = Réservoir souterrain (mm)
%
% FOLLOWING
%  - Boyle D.P. (2000) Multicriteria calibration of hydrological models, 
%    Ph.D. dissertation, Dep. of Hydrol. and Water Resour., Univ. of Arizona,
%    Tucson, USA
%  - Wagener, T., Boyle, D.P., Lees, M.J., Wheater, H.S., Gupta, H.V., 
%    Sorooshian, S., 2001. A framework for development and application of 
%    hydrological models. Hydrol. Earth Syst. Sci. 5, 13–26.
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

%% Modeling
% Get parameters
x = param.x ;
S = param.S ;
R1 = param.R1 ;
R2 = param.R2 ;
R3 = param.R3 ;
T = param.T ;
DL = param.DL ;
HY = param.HY ;

%%%
%%% PRODUCTION AND ROUTING
%%%
XS=S;
Ctprev = real(x(1)*(1-(1-(x(2)+1)*S/x(1))^(1/(x(2)+1))));
Ut1 = max(P-x(1)+Ctprev,0);
Pn = P-Ut1;
Dum = min(1,(Ctprev+Pn)/x(1));
S = x(1)/(x(2)+1)*(1-(1-Dum)^(x(2)+1));
Ut2 = max(0,Pn-(S-XS));
S = max(0,S-E);

Uq = x(3)*Ut2+Ut1;
Us = (1-x(3))*Ut2;
T = T+Us;
Qt = T/x(5)/x(6);
T = T-Qt;

R1 = R1+Uq;
Q1 = R1/x(6);
R1 = R1-Q1;

R2 = R2+Q1;
Q2 = R2/x(6);
R2 = R2-Q2;

R3 = R3+Q2;
Q3 = R3/x(6);
R3 = R3-Q3;

% Total discharge
%
HY = [HY(2:end); 0] + DL * (Qt+Q3) ;
Qsim = max( [0; HY(1)] ) ;

param.S = S;
param.R1 = R1;
param.R2 = R2;
param.R3 = R3;
param.T = T;
param.HY = HY ;

if nargout>2
    inter = [ S R1 R2 R3 T Qt Q1 Q2 Q3 ];
    
    % Flow components
    qsf = 0 ;
    qs1 = 0 ; qs2 = 0 ; qs = qs1+qs2 ;
    qrs1 = 0 ; qrs2 = 0 ; qrs = qrs1+qrs2 ;
    qss1 = 0 ; qss2 = 0 ; qss = qss1+qss2 ;
    qrss1 = Q3 ; qrss2 = 0 ; qrss = qrss1+qrss2 ;
    qn = Qt ;
    qr = 0 ;
    
    interq = [ qsf qs1 qs2 qs qrs1 qrs2 qrs qss1 qss2 qss qrss1 qrss2 qrss qn qr ];
    
    % Outputs
    varargout{1}=inter;
    varargout{2}=interq;
end