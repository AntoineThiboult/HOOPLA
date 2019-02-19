function [Qsim, param, varargout] = HydroMod15( P, E, param )
%
% [Qsim, param, varargout] = HydroMod15( P, E, param )
%
% SIMHYD hydrological model (modified version)
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
% inter = SIMHYD's internal values
% param -->
% 	.x(1) = Capacité du réservoir d'interception (mm)
% 	.x(2) = Capacité du réservoir de sol (mm)
% 	.x(3) = Constante de vidange du réservoir souterrain
% 	.x(4) = Délai
% 	.x(5) = Constante de vidange du réservoir de routage principal
% 	.x(6) = Constante d'écoulement hypodermique
% 	.x(7) = Constante de recharge de la nappe
% 	.x(8) = Infiltration maximale (mm)
%
% 	.S = Réservoir de sol
% 	.R = Réservoir souterrain
% 	.T = Réservoir de routage
%
% FOLLOWING
% Chiew, F.H.S., Peel, M.C., Western, A.W., 2002. Application and testing 
% of the simple rainfall-runoff model SIMHYD. pp. 335–367. (University of 
% Melbourne, Australia)
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

%%%
%%% PRODUCTION AND ROUTING
%%%

CAP1=x(1);
CAP=min(CAP1,E);
CAP=min(P,CAP);

if P > CAP
    EXC=P-CAP;
    E1=CAP;
else
    EXC=0;
    E1=P;
end

COEF=x(8);
SQ=2;
RINF=COEF*exp(-SQ*S/x(2));

if EXC > RINF
    SRUN=EXC-RINF;
    FILT=RINF;
else
    SRUN=0;
    FILT=EXC;
end

SINT=S/x(2)*FILT/x(6);
REC=max(0,S/x(2)*(FILT-SINT)/x(7));
S=S+FILT-SINT-REC;
EX2=0;

if S > x(2)
    EX2=S-x(2);
    S=x(2);
end

CAP2=10;
ET=min(E-E1,CAP2*S/x(2));
S=max(0,S-ET);

R=R+EX2+REC;
Qr=R/x(3)/x(5);
R=R-Qr;

T=T+SINT+SRUN+Qr;
Qt=T/x(5);
T=T-Qt;

% Total discharge
%
HY = [HY(2:end); 0] + DL * (Qt) ;
Qsim = max( [0; HY(1)] ) ;

param.S = S;
param.R = R;
param.T = T;
param.HY = HY ;

if nargout>2
    inter = [ S R T Qt Qr ];
    
    % Flow components
    qsf = 0 ;
    qs1 = 0 ; qs2 = 0 ; qs = qs1+qs2 ;
    qrs1 = 0 ; qrs2 = 0 ; qrs = qrs1+qrs2 ;
    qss1 = 0 ; qss2 = 0 ; qss = qss1+qss2 ;
    qrss1 = 0 ; qrss2 = 0 ; qrss = qrss1+qrss2 ;
    qn = 0 ;
    qr = Qt ;
    
    interq = [ qsf qs1 qs2 qs qrs1 qrs2 qrs qss1 qss2 qss qrss1 qrss2 qrss qn qr ];
    
    % Outputs
    varargout{1}=inter;
    varargout{2}=interq;
end