function [Qsim, param, varargout] = HydroMod1( P, E, param )
%
% [Qsim, param, varargout] = HydroMod1( P, E, param )
%
% BUCKET hydrological model
%
% INPUTS (time series of daily observations [n,1])
% P       = mean areal rainfall (mm)
% E       = mean areal evapotranspiration (mm)
% Q       = stream flow (mm)
% x       = the six model parameters (see "param" below) - [6,1]
%
% OUTPUTS
% Qs      = simulated stream flow (mm)
% perf    = model performances
% inter   = BUCKET's internal values
% param -->
%   .x(1) = capacité du réservoir sol
%   .x(2) = constante de dissociation du débordement du réservoir sol
%   .x(3) = constante de vidange du réservoir de routage (R)
%   .x(4) = délai
%   .x(5) = coefficient de partition de la pluie
%   .x(6) = constante de vidange du réservoir de routage (R,T)
%
%   .S    = Réservoir de sol
%   .R    = Réservoir de la couche racinaire (sous-sol)
%   .T    = Réservoir de routage direct
%
% FOLLOWING
%  - Thornthwaite, C.W., Mather, J.R., 1955. The water balance. Report. 
%    (Drexel Institute of Climatology. United States)
%  - Perrin, C. (2000). Vers une amélioration d'un modèle global pluie-débit, 
%    PhD Thesis, Appendix 1, p. 313-316. Retrieved from 
%    https://tel.archives-ouvertes.fr/tel-00006216
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

%% Modeling
% Get parameters
%
x   = param.x ;
S   = param.S ;
R   = param.R ;
T   = param.T ;
DL = param.DL ;
HY = param.HY ;

%%%
%%% PRODUCTION PART
%%%

% Soil moisture accounting (S)
%
Ps = (1-x(5))*P;
Pr = P-Ps;

if Ps >= E
    S = S+Ps-E;
    Is = max(0,S-x(1));
    S = S-Is;
else
    S = S*exp((Ps-E)/x(1));
    Is = 0;
end

%%%
%%% ROUTING PART
%%%

% Slow routing (R)
%
R = R+Is*(1-x(2));
Qr = R/(x(3)*x(6));
R = R-Qr;

% Fast routing (T)
%
T = T+Pr+Is*x(2);
Qt = T/x(6);
T = T-Qt;

% Calcul du débit total
%
HY = [HY(2:end); 0] + DL * (Qt+Qr) ;
Qsim = max( [0; HY(1)] ) ;


% Data
param.S  = S ;
param.R  = R ;
param.T  = T ;
param.HY = HY ;

if nargout>2
    inter = [ S R T Ps Pr Is Qr Qt ] ;
    
    % Flow components
    qsf = 0 ;
    qs1 = 0 ; qs2 = 0 ; qs = qs1+qs2 ;
    qrs1 = Qt ; qrs2 = 0 ; qrs = qrs1+qrs2 ;
    qss1 = 0 ; qss2 = 0 ; qss = qss1+qss2 ;
    qrss1 = 0 ; qrss2 = 0 ; qrss = qrss1+qrss2 ;
    qn = 0 ;
    qr = Qr ;
    
    interq = [ qsf qs1 qs2 qs qrs1 qrs2 qrs qss1 qss2 qss qrss1 qrss2 qrss qn qr ];
    
    % Outputs
    varargout{1}=inter;
    varargout{2}=interq;
end