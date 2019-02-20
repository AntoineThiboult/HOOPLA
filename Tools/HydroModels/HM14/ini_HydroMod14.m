function [Result, Param] = ini_HydroMod14(Switches, Date, x)
%
% [Result, Param] = ini_HydroMod14(Switches, Date, x)
%
% HydroMod14 initialization, for details see associated HydroMod14.m file
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

% Parameters
Param.x = x ;

% Routing delay consideration
drtfc = ceil( x(9) );
k = (0:drtfc)';
Param.DL = zeros( size(k) );
Param.DL(end-1) = 1/(x(9)-k(end-1)+1);
Param.DL(end) = 1-1/(x(9)-k(end-1)+1);
Param.HY  = zeros( size(Param.DL) );

% Fixed parameters
Param.XF1 = 3;
Param.XF2 = 30;

% Initialization of the reservoir states
%
Param.S = 3;
Param.T = 10;
Param.R = 100;
Param.L = 0;
Param.M = 0;

% Initialize HydroMod14 for all time steps
%
lP    = length( Date ) ;
Result.Qs    = zeros( lP,1 ) ;

if Switches.exportLight.on == 0
    % Flow components
    Result.inter = zeros( lP, 11 ) ;
    Result.interq = zeros( lP, 15 ) ;
end