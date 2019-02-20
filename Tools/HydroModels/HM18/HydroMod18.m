function [Qsim, param, varargout] = HydroMod18( P, E, param )
% 
% [Qsim, param, varargout] = HydroMod18( P, E, param )
%
% TOPMODEL hydrological model (modified version)
%
% INPUTS
% P 	= mean areal rainfall (mm)
% E 	= mean areal evapotranspiration (mm)
% param	= the ten model parameters (see "param" below) - [7,1]
%
% OUTPUTS
% Qsim    = simulated stream flow (mm)
% inter   = HydroMod18's internal values (varargout 1)
% interq  = HydroMod18's internal flow components (varargout 2)
% param -->
% 	.x(1) = Quadratic routing reservoir capacity
% 	.x(2) = Exponential emptying paramter
% 	.x(3) = Interception reservoir capacity
% 	.x(4) = Delay
% 	.x(5) = Topographic coefficient distribution parameter
% 	.x(6) = Topographic index parameter
% 	.x(7) = PET parameter
%
% 	.S = Interception reservoir state (mm)
% 	.T = Ground reservoir state (mm)
% 	.R = Quadratic routing reservoir state (mm)
%
% FOLLOWING
%  - Beven, K.J., Kirkby, M.J., Schofield, N., Tagg, A.F., 1984. Testing a 
%    physically-based flood forecasting model (TOPMODEL) for 3 UK catchments. 
%    J. Hydrol. 69, 119�143. (Institute of Environmental and Biological
%    Sciences, University of Lancaster, UK)
%  - Perrin, C. (2000). Vers une am�lioration d'un mod�le global pluie-d�bit, 
%    PhD Thesis, Appendix 1, p. 452-458. Retrieved from
%    https://tel.archives-ouvertes.fr/tel-00006216
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

%% Modeling
% Get parameters
%
x = param.x ;
S = param.S ;
T = param.T ;
R = param.R ;
DL = param.DL ;
HY = param.HY ;
T=T-1e4; % ajoutAT for EnKF issues. 

%%% PRODUCTION PART
%
% Inteception storage (S)
%
S = S+P;
Es = min(S,E);
S = S-Es;
Er = E-Es;
Pr = max(0,S-x(3));
S = S-Pr;
Ps = Pr / (1+exp(x(6)-(T/x(5))));

% Groundwater storage (T)
%
T = T+Pr-Ps;
Et = Er / (1+exp(x(7)-(T/x(5))));
T = T+Et;

%%% ROUTING
%
% Routing storage (R)
%
R = R+Ps;
Qr = (R^2)/(R+x(1));
R = R-Qr;

% Baseflow
Qt = x(2)*exp(T/x(2));
T = T-Qt;

% Total discharge
HY = [HY(2:end); 0] + DL * (Qt+Qr) ;
Qsim = max( [0; HY(1)] ) ;

T=T+1e4; % ajoutAT for EnKF issues. 
param.T = T;
param.S = S;
param.R = R;
param.HY = HY ;

if nargout>2
    inter = [ S T R Pr Ps Qt Qr ];
    
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