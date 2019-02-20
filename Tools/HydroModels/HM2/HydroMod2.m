function [Qsim, param, varargout] = HydroMod2( P, E, param )
%
% [Qsim, param, varargout] = HydroMod2( P, E, param )
%
% CEQUEAU hydrological model (modified version)
%
% INPUTS 
% P 	= mean areal rainfall (mm)
% E 	= mean areal evapotranspiration (mm)
% param = the nine model parameters (see "param" below) - [9,1]
%
% OUTPUTS
% Qsim    = simulated stream flow (mm)
% inter   = HydroMod2's internal values (varargout 1)
% interq  = HydroMod2's internal flow components (varargout 2)
% param -->
% 	.x(1) = infiltratin threshold (mm)
% 	.x(2) = 1rst reservoir emptying threshold (mm)
% 	.x(3) = emptying infiltration constant
% 	.x(4) = upper lateral emptying constant of soil reservoir
% 	.x(5) = maximum capacity of ground reservoir (mm)
% 	.x(6) = delay
% 	.x(7) = ground reservoir emptying threshold (mm)
% 	.x(8) = lower lateral emptying constant of soil reservoir
% 	.x(9) = lower lateral emptying constant of ground reservoir
%
% 	.S = Surface reservoir state (mm)
% 	.T = Ground reservoir state (mm)
%
% FOLLOWING
%  - Girard, G., Morin, G., Charbonneau, R., 1972. Modèle précipitations-débits
%    à discrétisation spatiale. Cahiers ORSTOM, Série Hydrologie 9, 35–52.
%    (INRS Eau, Québec, Canada)
%  - Perrin, C. (2000). Vers une amélioration d'un modèle global pluie-débit,
%    PhD Thesis, Appendix 1, p. 322-326. Retrieved from
%    https://tel.archives-ouvertes.fr/tel-00006216
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
