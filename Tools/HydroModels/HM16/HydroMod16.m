function [Qsim, param, varargout] = HydroMod16( P, E, param )
%
% [Qsim, param, varargout] = HydroMod16( P, E, param )
%
% SMAR hydrological model (modified version)
%
% INPUTS (time series of daily observations [n,1])
% P 	= mean areal rainfall (mm)
% E 	= mean areal evapotranspiration (mm)
% Q 	= stream flow (mm)
% x 	= the eight model parameters (see "param" below) - [8,1]
%
% OUTPUTS
% Qs    = simulated stream flow (mm)
% perf  = model performances
% inter = SMAR's internal values
% param -->
% 	.x(1) = paramètre d’écoulement direct
% 	.x(2) = paramètre d’infiltration Yc
% 	.x(3) = coefficient de réduction de l’ETP C
% 	.x(4) = capacité du réservoir quadratique
% 	.x(5) = constante de vidange du réservoir de routage linéaire
% 	.x(6) = délai
% 	.x(7) = paramètre de correction de l’ETP
% 	.x(8) = coefficient de partage G
%
% 	.S = Réservoir de sol multicouches (N couches, de hauteur XF1)(mm)
% 	.L = Réservoir de routage linéaire (mm)
% 	.T = Réservoir de routage quadratique (mm)
%
% FOLLOWING
% O’Connell et al. (1970)
% Department of Engineering Hydrology, University College, Galway, Ireland
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

%% Modeling
% Get parameters
%
x = param.x ;
S = param.S ;
L = param.L ;
T = param.T ;
XF1 = param.XF1 ;
XF2 = param.XF2 ;
s1 = param.s1;
s2 = param.s2;
s3 = param.s3;
s4 = param.s4;
s5 = param.s5;
s6 = param.s6;
s7 = param.s7;
s8 = param.s8;
s9 = param.s9;
s10 = param.s10;
s11 = param.s11;
s12 = param.s12;
s13 = param.s13;
s14 = param.s14;
s15 = param.s15;
s16 = param.s16;
DL = param.DL ;
HY = param.HY ;

% Net inputs
%
E = x(7)*E;
Pn = max(0,P-E);
En = max(0,E-P);

% Direct runoff and infiltration
%
Pr1 = x(1)*(S/(5*XF1))*Pn;
Fr = XF2*exp((-x(2))*(S/(5*XF1)));
Ps = min(Fr, Pn-Pr1);
Pr2 = Pn-Pr1-Ps;

% Multi-layers soil moisture accounting (S)
%
s1 = s1+Ps;
Ps = max(0,s1-XF1);
s1 = s1-Ps;
Es1 = min(s1,(x(3)^1)*En);
s1 = s1-Es1;
En = En-Es1;

s2 = s2+Ps;
Ps = max(0,s2-XF1);
s2 = s2-Ps;
Es2 = min(s2,(x(3)^2)*En);
s2 = s2-Es2;
En = En-Es2;

s3 = s3+Ps;
Ps = max(0,s3-XF1);
s3 = s3-Ps;
Es3 = min(s3,(x(3)^3)*En);
s3 = s3-Es3;
En = En-Es3;

s4 = s4+Ps;
Ps = max(0,s4-XF1);
s4 = s4-Ps;
Es4 = min(s4,(x(3)^4)*En);
s4 = s4-Es4;
En = En-Es4;

s5 = s5+Ps;
Ps = max(0,s5-XF1);
s5 = s5-Ps;
Es5 = min(s5,(x(3)^5)*En);
s5 = s5-Es5;
En = En-Es5;

s6 = s6+Ps;
Ps = max(0,s6-XF1);
s6 = s6-Ps;
Es6 = min(s6,(x(3)^6)*En);
s6 = s6-Es6;
En = En-Es6;

s7 = s7+Ps;
Ps = max(0,s7-XF1);
s7 = s7-Ps;
Es7 = min(s7,(x(3)^7)*En);
s7 = s7-Es7;
En = En-Es7;

s8 = s8+Ps;
Ps = max(0,s8-XF1);
s8 = s8-Ps;
Es8 = min(s8,(x(3)^8)*En);
s8 = s8-Es8;
En = En-Es8;

s9 = s9+Ps;
Ps = max(0,s9-XF1);
s9 = s9-Ps;
Es9 = min(s9,(x(3)^9)*En);
s9 = s9-Es9;
En = En-Es9;

s10 = s10+Ps;
Ps = max(0,s10-XF1);
s10 = s10-Ps;
Es10 = min(s10,(x(3)^10)*En);
s10 = s10-Es10;
En = En-Es10;

s11 = s11+Ps;
Ps = max(0,s11-XF1);
s11 = s11-Ps;
Es11 = min(s11,(x(3)^11)*En);
s11 = s11-Es11;
En = En-Es11;

s12 = s12+Ps;
Ps = max(0,s12-XF1);
s12 = s12-Ps;
Es12 = min(s12,(x(3)^12)*En);
s12 = s12-Es12;
En = En-Es12;

s13 = s13+Ps;
Ps = max(0,s13-XF1);
s13 = s13-Ps;
Es13 = min(s13,(x(3)^13)*En);
s13 = s13-Es13;
En = En-Es13;

s14 = s14+Ps;
Ps = max(0,s14-XF1);
s14 = s14-Ps;
Es14 = min(s14,(x(3)^14)*En);
s14 = s14-Es14;
En = En-Es14;

s15 = s15+Ps;
Ps = max(0,s15-XF1);
s15 = s15-Ps;
Es15 = min(s15,(x(3)^15)*En);
s15 = s15-Es15;
En = En-Es15;

s16 = s16+Ps;
Ps = max(0,s16-XF1);
s16 = s16-Ps;
Es16 = min(s16,(x(3)^16)*En);
s16 = s16-Es16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = s1+s2+s3+s4+s5;

Is = Ps;

% Linear routing storage (L)
%
L = L+(1-x(8))*Is;
Ql = L/x(5);
L = L-Ql;

% Quadratic routing storage (T)
%
T = T+(x(8)*Is)+Pr2;
Qt = (T^2)/(T+x(4));
T = T-Qt;

% Total discharge
%
HY = [HY(2:end); 0] + DL * (Ql+Qt+Pr1) ;
Qsim = max( [0; HY(1)] ) ;

param.S = S;
param.L = L;
param.T = T;
param.s1 = s1;
param.s2 = s2;
param.s3 = s3;
param.s4 = s4;
param.s5 = s5;
param.s6 = s6;
param.s7 = s7;
param.s8 = s8;
param.s9 = s9;
param.s10 = s10;
param.s11 = s11;
param.s12 = s12;
param.s13 = s13;
param.s14 = s14;
param.s15 = s15;
param.s16 = s16;
param.HY = HY ;

if nargout>2
    inter = [ S L T Ps Pr1 Pr2 Ql Qt ];
    
    % Flow components
    qsf = Pr1 ;
    qs1 = 0 ; qs2 = 0 ; qs = qs1+qs2 ;
    qrs1 = 0 ; qrs2 = 0 ; qrs = qrs1+qrs2 ;
    qss1 = 0 ; qss2 = 0 ; qss = qss1+qss2 ;
    qrss1 = 0 ; qrss2 = 0 ; qrss = qrss1+qrss2 ;
    qn = Ql ;
    qr = Qt ;
    
    interq = [ qsf qs1 qs2 qs qrs1 qrs2 qrs qss1 qss2 qss qrss1 qrss2 qrss qn qr ];
    
    % Outputs
    varargout{1}=inter;
    varargout{2}=interq;
end