function [Result, Param] = ini_HydroMod8(Switches, Date, x)
%
% [Result, Param] = ini_HydroMod8(Switches, Date, x)
%
% IHACRES initialization, for details see associated IHACRES.m file
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

% Parameters
Param.x = x ;

% Routing delay consideration
drtfc = ceil( x(5) );
k = (0:drtfc)';
Param.DL = zeros( size(k) );
Param.DL(end-1) = 1/(x(5)-k(end-1)+1);
Param.DL(end) = 1-1/(x(5)-k(end-1)+1);
Param.HY  = zeros( size(Param.DL) );

% Initialization of the reservoir states
%
Param.S = 0.5;
Param.T = 50;
Param.R = 5;

% Apply IHACRES for all time steps
%
lP    = length( Date ) ;
Result.Qs    = zeros( lP,1 ) ;

if Switches.exportLight.on == 0
    % Flow components
    Result.inter = zeros( lP, 8 ) ;
    Result.interq = zeros( lP, 15 ) ;
end