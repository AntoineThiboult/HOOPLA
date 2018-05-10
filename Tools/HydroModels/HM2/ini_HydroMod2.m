function [Result, Param] = ini_HydroMod2(Switches, Date, x)
%
% [Result, Param] = ini_HydroMod2(Switches, Date, x)
%
% CEQUEAU initialization, for details see associated CEQUEAU.m file
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
Param.S = 500;
Param.T = x(5)*0.2;

% Apply CEQUEAU for all time steps
%
lP    = length( Date ) ;
Result.Qs    = zeros( lP,1 ) ;

if Switches.exportLight.on == 0
    % Flow components
    Result.inter = zeros( lP, 13 ) ;
    Result.interq = zeros( lP, 15 ) ;
end