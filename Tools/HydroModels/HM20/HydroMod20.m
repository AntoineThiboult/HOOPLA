function [Qsim, param, varargout] = HydroMod20( P, E, param )
%
% [Qsim, param, varargout] = HydroMod20( P, E, param )
%
% XINANJIANG hydrological model (modified version)
%
% INPUTS
% P 	= mean areal rainfall (mm)
% E 	= mean areal evapotranspiration (mm)
% param	= the eight model parameters (see "param" below) - [8,1]
%
% OUTPUTS
% Qsim    = simulated stream flow (mm)
% inter   = HydroMod20's internal values (varargout 1)
% interq  = HydroMod20's internal flow components (varargout 2)
% param -->
% 	.x(1) = Flow partioning coefficient
% 	.x(2) = Fast emptying reservoir coefficient
% 	.x(3) = Slow emptying reservoir coefficient
% 	.x(4) = Water-soil reservoir capacity
% 	.x(5) = Soil reservoir capacity
% 	.x(6) = Delay
% 	.x(7) = Water-soil reservoir emptying coefficient
% 	.x(8) = Free-water reservoir exponent
%
% 	.S = Soil reservoir state (mm)
% 	.R = Water-soil reservoir state (mm)
% 	.T = Fast routing reservoir state (mm)
% 	.M = Slow routing reservoir state (mm)
%
% FOLLOWING
%  - Zhao et al. (1980), Zhao and Liu (1995), The Xinanjiang model. IAHS
%    Publications 129, 351356. (Ho-hai University, Nanjing, China)
%  - Perrin, C. (2000). Vers une amélioration d'un modèle global pluie-débit, 
%    PhD Thesis, Appendix 1, p. 463-467. Retrieved from 
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
M = param.M ;
XF1 = 1/4;
DL = param.DL ;
HY = param.HY ;

%%% PRODUCTION PART

% Net inputs
%
Pn = max(0,P-E);
En = max(0,E-P);

% If Evaporation
%
%#ok<*NASGU>
if En > 0
    if S/x(5) >= 0.9
        Es = min(S,En);
    elseif S/x(5) < 0.09
        Es = min(S,En*0.1);
    else
        Es = min(S,(En*S)/(0.9*x(5)));
    end
    S=S-Es;
    Qs0 = 0; 
    Ir = 0;
else
    En = 0;
    Ir = 0;
end

% If Net rainfall
%
if Pn > 0
    Fs = ( ( max(0,1-(S/x(5))) )^(1/(1+XF1)) )-(Pn/((1+XF1)*x(5)));
    Fs = (max(Fs,0))^(1+XF1);
    Ps = max(0,x(5)-S-Fs*x(5));
    S = min(x(5),S+Ps);
    Pr = max(0,Pn-Ps);
    Fr = ( (max(0,1-(R/x(4))))^(1/(1+x(8))) )-(Pr/((1+x(8))*x(4)));
    Fr = (max(Fr,0))^(1+x(8));
    Pr2 = max(0,x(4)-R-Fr*x(4));
    R = min(x(4),R+Pr2);
    Qs0 = max(0,Pr-Pr2);
    Ir = R/x(7);
else
    Pn = 0;
    Ir = 0;
    Pr = 0;
    Ps = 0;
    Qs0 = 0;
end

%%% UPDATE and ROUTING
%
R = R-Ir;
T = T+Ir*x(1);
Qt = T/x(2);
T = T-Qt;
M = M+Ir*(1-x(1));
Qm = M/(x(2)*x(3));
M = M-Qm;

% Total discharge
HY = [HY(2:end); 0] + DL * (Qs0 + Qm + Qt) ;
Qsim = max( [0; HY(1)] ) ;

param.S = S;
param.R = R;
param.T = T;
param.M = M;
param.HY = HY;

if nargout>2
    inter = [ S R T M Pn Pr Ps Qs0 Qt Qm ];
    
    % Flow components
    qsf = Qs0 ;
    qs1 = 0 ; qs2 = 0 ; qs = qs1+qs2 ;
    qrs1 = 0 ; qrs2 = 0 ; qrs = qrs1+qrs2 ;
    qss1 = 0 ; qss2 = 0 ; qss = qss1+qss2 ;
    qrss1 = 0 ; qrss2 = 0 ; qrss = qrss1+qrss2 ;
    qn = Qm ;
    qr = Qt ;
    
    interq = [ qsf qs1 qs2 qs qrs1 qrs2 qrs qss1 qss2 qss qrss1 qrss2 qrss qn qr ];
    
    % Outputs
    varargout{1}=inter;
    varargout{2}=interq;
end