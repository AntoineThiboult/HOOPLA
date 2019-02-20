function [Result, Param] = ini_HydroMod6(Switches, Date, x)
%
% [Result, Param] = ini_HydroMod6(Switches, Date, x)
%
% HydroMod6 initialization, for details see associated HydroMod6.m file
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

% Parameters
Param.x = x ;

% Initialization of the reservoir states
%
Param.S = x(1);
Param.R = 1;
Param.T = 10;

% UH initialization (triangular weighting function)
%
t1=1:floor(x(6)/2); % Time of rising 
t2=t1(end)+1:ceil(x(6)); % Time of recession
UH= [t1-0.5, x(6)+0.5-t2]'; % Triangular hydrograph
Param.UH = UH./sum(UH); % Normalization
Param.H  = zeros( size(Param.UH) ) ;

% Initialize HydroMod6 for all time steps
%
lP    = length( Date ) ;
Result.Qs    = zeros( lP,1 ) ;

if Switches.exportLight.on == 0
    % Flow components
    Result.inter = zeros( lP, 7 ) ;
    Result.interq = zeros( lP, 15 ) ;
end