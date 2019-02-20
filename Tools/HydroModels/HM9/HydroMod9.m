function [Qsim, param, varargout] = HydroMod9( P, E, param )
%
% [Qsim, param, varargout] = HydroMod9( P, E, param )
%
% MARTINE hydrological model (modified version)
%
% INPUTS
% P       = mean areal rainfall (mm)
% E       = mean areal evapotranspiration (mm)
% param   = the seven model parameters (see "param" below) - [7,1]
%
% OUTPUTS
% Qsim    = simulated stream flow (mm)
% inter   = HydroMod9's internal values (varargout 1)
% interq  = HydroMod9's internal flow components (varargout 2)
% param -->
%   .x(1) = Surface reservoir capacity
%   .x(2) = Intermediate reservoir capacity
%   .x(3) = Quadratic routing reservoir capacity
%   .x(4) = Emptying constant of ground reservoir
%   .x(5) = Distribution coefficient
%   .x(6) = Delay
%   .x(7) = Intermediate reservoir emptying constant
%   .S    = Surface reservoir state
%   .T    = Intermediate reservoir state
%   .L    = Ground reservoir state
%   .R    = Quadratic routing reservoir state
%
% FOLLOWING
%  - Mazenc, B., Sanchez, M., Thiery, D., 1984. Analyse de la influence de 
%    la physiographie d’un bassin versant sur les paramètres d’un modèle
%    hydrologique global et sur les débits caractéristiques à la exutoire. 
%    J. Hydrol. 69, 97–188. (Bureau de Recherches Géologiques et Minières, 
%    Orléans, France)
%  - Perrin, C. (2000). Vers une amélioration d'un modèle global pluie-débit, 
%    PhD Thesis, Appendix 1, p. 387-390. Retrieved from 
%    https://tel.archives-ouvertes.fr/tel-00006216
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

%% Modeling
% Get parameters
%
x   = param.x ;
S   = param.S ;
T   = param.T ;
L   = param.L ;
R   = param.R ;
DL = param.DL ;
HY = param.HY ;

%%%
%%% PRODUCTION AND ROUTING
%%%

% Soil storage (S)
%
S = S+P;
Pr = max(0,S-x(1));
S = S-Pr;
Es = min(S,E);
S = S-Es;
Et = E-Es;

% Direct routing storage (R)
%
R = R+x(5)*Pr;
Qr = (R^2)/(R+x(3));
R = R-Qr;

% Lower storage (T)
%
T = max(0,T-Et);
T = T+(1-x(5))*Pr;
Qt1 = max(0,T/x(7));
T = T-Qt1;
Qt2 = max(0,T-x(2));
T = T-Qt2;

% Groundwater storage (L)
%
L = L+Qt1+Qt2;
Ql = L/x(4);
L = L-Ql;

% Total discharge
%
HY = [HY(2:end); 0] + DL * (Ql+Qr) ;
Qsim = max( [0; HY(1)] ) ;

% Data
param.S  = S ;
param.T  = T ;
param.L  = L ;
param.R  = R ;
param.HY = HY ;

if nargout>2
    inter = [ S T L R Pr Es Et Qt1 Qt2 Ql Qr ] ;
    
    % Flow components
    qsf = 0 ;
    qs1 = 0 ; qs2 = 0 ; qs = qs1+qs2 ;
    qrs1 = 0 ; qrs2 = 0 ; qrs = qrs1+qrs2 ;
    qss1 = 0 ; qss2 = 0 ; qss = qss1+qss2 ;
    qrss1 = 0 ; qrss2 = 0 ; qrss = qrss1+qrss2 ;
    qn = Ql ;
    qr = Qr ;
    
    interq = [ qsf qs1 qs2 qs qrs1 qrs2 qrs qss1 qss2 qss qrss1 qrss2 qrss qn qr ];
    
    % Outputs
    varargout{1}=inter;
    varargout{2}=interq;
end