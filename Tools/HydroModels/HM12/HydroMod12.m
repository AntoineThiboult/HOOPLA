function [Qsim, param, varargout] = HydroMod12( P, E, param )
%
% [Qsim, param, varargout] = HydroMod12( P, E, param )
%
% NAM10 hydrological model (modified version)
%
% INPUTS
% P 	= mean areal rainfall (mm)
% E 	= mean areal evapotranspiration (mm)
% param	= the ten model parameters (see "param" below) - [10,1]
%
% OUTPUTS
% Qsim    = simulated stream flow (mm)
% inter   = HydroMod12's internal values (varargout 1)
% interq  = HydroMod12's internal flow components (varargout 2)
% param -->
% 	.x(1) = Emptying threshold of the ground reservoir
% 	.x(2) = Emptying constant of routing reservoirs
% 	.x(3) = Sub-surface flow constant
% 	.x(4) = Delay
% 	.x(5) = Percolation constant
% 	.x(6) = Emptying constant of the ground reservoir
% 	.x(7) = Maximum capacity of the surface reservoir
% 	.x(8) = Surface flow constant
% 	.x(9) = Maximum capacity of soil reservoir
% 	.x(10) = Capillary rise parameter
%
% 	.U = Surface reservoir state (mm)
% 	.L = Soil reservoir state (mm)
% 	.CK1 = Overlandflow 1 reservoir state (mm)
% 	.CK2 = Overlandflow 2 reservoir state (mm)
% 	.CK1b = Interflow 1 reservoir state(mm)
% 	.CK2b = Interflow 2 reservoir state(mm)
% 	.GW = Ground water reservoir state (mm)
%
% FOLLOWING
%  - Nielsen, S.A., Hansen, E., 1973. Numerical simulation of the 
%    rainfall-runoff process on a daily basis. Nord. Hydrol. 4, 171–190.
%    (Institute of Hydrodynamics and Hydraulic Engineering, Technical 
%    University of Denmark, and The Danish Hydraulic Institute of Denmark.)
%  - Perrin, C. (2000). Vers une amélioration d'un modèle global pluie-débit, 
%    PhD Thesis, Appendix 1, p. 411-415. Retrieved from 
%    https://tel.archives-ouvertes.fr/tel-00006216
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

%% Modeling
% Get parameters
%
x = param.x ;
U = param.U ;
L = param.L ;
CK1 = param.CK1 ;
CK1b = param.CK1b ;
CK2 = param.CK2 ;
CK2b = param.CK2b ;
GW = param.GW ;
DL = param.DL ;
HY = param.HY ;

%%% PRODUCTION PART

U=U+P;
QIF=min(U,L/x(7)*U/x(3));
U=U-QIF;
CK1=CK1+QIF;
B21=CK1/x(2);
CK1=CK1-B21;
CK2=CK2+B21;
B2=CK2/x(2);
CK2=CK2-B2;

U=U-E;

if U >= 0 && U <= x(9)
    E1 = 0;
    PN = 0;
elseif U < 0
    E1 = -U;
    U = 0;
    PN = 0;
elseif U > x(9)
    E1 = 0;
    PN = U-x(9);
    U = x(9);
end

QOF=0;
DL0=0;
G=0;
if PN > 0
    QOF=PN*L/x(7)/x(8);
    if x(5) == 1
        G = 0;
    end
    if x(5) ~= 1
        if L/x(7) > x(5)
            G=(PN-QOF)*(L/x(7)-x(5))/(1-x(5));
        else
            G = 0;
        end
    end
    DL0=PN-QOF-G;
    if DL0 > x(7)
        DL1=DL0-(x(7)-L);
        G=G+DL1;
    end
end

CK1b=CK1b+QOF;
B12=CK1b/x(2);
CK1b=CK1b-B12;

CK2b=CK2b+B12;
B1=CK2b/x(2);
CK2b=CK2b-B1;

L=L+DL0;
L=max(0.,L-E1*L/x(7));

GW=GW-G;
if GW <= x(1)
    BF=(x(1)-GW)/x(6);
else
    BF=0;
end

GW=GW+BF;

if GW > 0
    BF1=0;
else
    BF1=-GW+0.1;
    GW=0.1;
end

if L > x(7)
    L = x(7);
end

CAFLU=((1-L/x(7))^0.5)*(x(10)/GW)^2;

if CAFLU > x(7)-L
    CAFLU = x(7)-L;
end

L=L+CAFLU;
GW=GW+CAFLU;

% Total discharge
%
HY = [HY(2:end); 0] + DL * (BF+BF1+B1+B2) ;
Qsim = max( [0; HY(1)] ) ;

param.U = U;
param.L = L;
param.CK1 = CK1;
param.CK1b = CK1b;
param.CK2 = CK2;
param.CK2b = CK2b;
param.GW = GW;
param.HY = HY ;

if nargout>2
    inter = [ U L CK1 CK2 CK1b CK2b GW ];
    
    % Flow components
    qsf = 0 ;
    qs1 = 0 ; qs2 = 0 ; qs = qs1+qs2 ;
    qrs1 = 0 ; qrs2 = B1 ; qrs = qrs1+qrs2 ;
    qss1 = 0 ; qss2 = 0 ; qss = qss1+qss2 ;
    qrss1 = 0 ; qrss2 = B2 ; qrss = qrss1+qrss2 ;
    qn = BF+BF1 ;
    qr = 0 ;
    
    interq = [ qsf qs1 qs2 qs qrs1 qrs2 qrs qss1 qss2 qss qrss1 qrss2 qrss qn qr ];
    
    % Outputs
    varargout{1}=inter;
    varargout{2}=interq;
end