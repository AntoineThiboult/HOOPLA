function [Qsim, param, varargout] = HydroMod10( P, E, param )
%
% [Qsim, param, varargout] = HydroMod10( P, E, param )
%
% MOHYSE hydrological model (modified version)
%
% INPUTS
% P       = mean areal rainfall (mm)
% E       = mean areal evapotranspiration (mm)
% param   = the seven model parameters (see "param" below) - [7,1]
%
% OUTPUTS
% Qsim    = simulated stream flow (mm)
% inter   = HydroMod10's internal values (varargout 1)
% interq  = HydroMod10's internal flow components (varargout 2)
% param -->
%   .x(1) = transpiration coefficient (Ctr)
%   .x(2) = maximum infiltration rate (Cinf)
%   .x(3) = Emptying coefficient of the aquifer vadose zone (Cva)
%   .x(4) = Emptying coefficient of the river vadose zone (Cv)
%   .x(5) = Emptying coefficient from aquifer to river coefficient (Ca)
%   .x(6) = Unit hydrograph shape parameter (alpha)
%   .x(7) = Unit hydrograph scale parameter (beta)
%   .S    = Soil reservoir state 
%   .R    = Ground reservoir state
%   .UH   = Unit hydrograph
%   .HU   = Hydrograph values (mm) - updated at each time step
%
% FOLLOWING
% Fortin, V., Turcotte, R., 2007. Le modèle hydrologique mohyse, note de
% cours pour SCA7420. Report, Département des sciences de la terre et de
% l’atmosphère, Université du Québec à Montreal.
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

%% Modeling
% Get parameters
%
x   = param.x ;
S   = param.S ;
R   = param.R ;
UH  = param.UH ;
HUQ1  = param.HUQ1 ;
HUQ2  = param.HUQ2 ;
HUQ3  = param.HUQ3 ;

%%%
%%% PRODUCTION AND ROUTING
%%%

ED = min(P,E);
TR = min(x(1)*S,E-ED);
if S >= x(2)
    I = 0;
else
    I = (P-ED)*(1-S/x(2));
end

% Calcul des flux

Q1 = P-ED-I;
Q2 = x(4)*S;
Q3 = x(5)*R;
qt = x(3)*S;

% Mise à jour des réservoirs

S = max(0,S+I-TR-qt-Q2);
R = R+qt-Q3;

% Mise à jour de l'hydrogramme
%
HUQ1 = [HUQ1(2:end); 0] + UH .* Q1 ;
HUQ2 = [HUQ2(2:end); 0] + UH .* Q2 ;
HUQ3 = [HUQ3(2:end); 0] + UH .* Q3 ;


% Calcul du débit total
%
Qsim = HUQ1(1)+HUQ2(1)+HUQ3(1);

param.S  = S ;
param.R  = R ;
param.HUQ1 = HUQ1 ;
param.HUQ2 = HUQ2 ;
param.HUQ3 = HUQ3 ;

if nargout>2
    inter = [ S R ED I TR Q1 Q2 Q3 ] ;
    
    % Flow components
    qsf = Q1 ;
    qs1 = Q2 ; qs2 = 0 ; qs = qs1+qs2 ;
    qrs1 = 0 ; qrs2 = 0 ; qrs = qrs1+qrs2 ;
    qss1 = 0 ; qss2 = 0 ; qss = qss1+qss2 ;
    qrss1 = 0 ; qrss2 = 0 ; qrss = qrss1+qrss2 ;
    qn = Q3;
    qr = 0 ;
    
    interq = [ qsf qs1 qs2 qs qrs1 qrs2 qrs qss1 qss2 qss qrss1 qrss2 qrss qn qr ];
    
    % Outputs
    varargout{1}=inter;
    varargout{2}=interq;
end