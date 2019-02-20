function [Result, Param] = ini_HydroMod9(Switches, Date, x)
%
% [Result, Param] = ini_HydroMod9(Switches, Date, x)
%
% HydroMod9 initialization, for details see associated HydroMod9.m file
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

% Parameters
Param.x = x ;

% Routing delay consideration
drtfc = ceil( x(6) );
k = (0:drtfc)';
Param.DL = zeros( size(k) );
Param.DL(end-1) = 1/(x(6)-k(end-1)+1);
Param.DL(end) = 1-1/(x(6)-k(end-1)+1);
Param.HY  = zeros( size(Param.DL) );

% Initialization of the reservoir states
%
Param.S = x(1) ;
Param.T = x(2)*0.5 ;
Param.L = 5 ;
Param.R = x(3)*0.1 ;

% Initialize HydroMod9 for all time steps
%
lP    = length( Date ) ;
Result.Qs    = zeros( lP,1 ) ;

if Switches.exportLight.on == 0
    % Flow components
    Result.inter = zeros( lP, 11 ) ;
    Result.interq = zeros( lP, 15 ) ;
end