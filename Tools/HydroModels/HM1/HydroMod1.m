function [Qsim, param, varargout] = HydroMod1( P, E, param )
%
% [Qsim, param, varargout] = HydroMod1( P, E, param )
%
% BUCKET hydrological model
%
% INPUTS
% P       = mean areal rainfall (mm)
% E       = mean areal evapotranspiration (mm)
% param   = the six model parameters (see "param" below) - [6,1]
%
% OUTPUTS
% Qsim    = simulated stream flow (mm)
% inter   = HydroMod1's internal values (varargout 1)
% interq  = HydroMod1's internal flow components (varargout 2)
% param -->
%   .x(1) = soil reservoir capacity.
%   .x(2) = soil reservoir overflow dissociation constant
%   .x(3) = routing reservoir emptying constant (reservoir R)
%   .x(4) = delay
%   .x(5) = rainfall partitioning coefficient
%   .x(6) = routing reservoir emptying constant (reservoir R and T)
%
%   .S    = Soil reservoir state
%   .R    = Root layer reservoir state
%   .T    = Direct routing reservoir state
%
% FOLLOWING
%  - Thornthwaite, C.W., Mather, J.R., 1955. The water balance. Report. 
%    (Drexel Institute of Climatology. United States)
%  - Perrin, C. (2000). Vers une am�lioration d'un mod�le global pluie-d�bit, 
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

% Calcul du d�bit total
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