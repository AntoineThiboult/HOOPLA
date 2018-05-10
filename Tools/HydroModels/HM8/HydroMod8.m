function [Qsim, param, varargout] = HydroMod8( P, E, param )
%
% [Qsim, param, varargout] = HydroMod8( P, E, param )
%
% IHACRES hydrological model (modified version)
%
% INPUTS (time series of daily observations [n,1])
% P 	= mean areal rainfall (mm)
% E 	= mean areal evapotranspiration (mm)
% Q 	= stream flow (mm)
% x 	= the seven model parameters (see "param" below) - [7,1]
%
% OUTPUTS
% Qs    = simulated stream flow (mm)
% perf  = model performances
% inter = IHACRES's internal values
% param -->
% 	.x(1) = Param�tre de for�age 1/C (calcul�e de sorte que le volume de pluie efficace calcul� soit �gal au volume des apports en eau observ�s sur la p�riode de calage)
% 	.x(2) = Param�tre de partage des �coulements
% 	.x(3) = Constante de vidange du r�servoir de routage rapide
% 	.x(4) = Constante de vidange du r�servoir de routage lent (X3.X4)
% 	.x(5) = D�lai
% 	.x(6) = Param�tre f (facteur de modulation de la temp�rature)
% 	.x(7) = Param�tre Tw (cste caract�ristique de l'ass�chement du bassin)
%
% 	.S = R�servoir de suivi de l'humidit� (mm)
% 	.T = R�servoir de routage lent (mm)
% 	.R = R�servoir de routage rapide (mm)
%
% FOLLOWING
% Littlewood et al. (1997)
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

%% Modeling
% Get parameters
%
x   = param.x ;
S   = param.S ;
T   = param.T ;
R   = param.R ;
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