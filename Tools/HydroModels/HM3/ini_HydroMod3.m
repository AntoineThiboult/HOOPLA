function [Result, Param] = ini_HydroMod3(Switches, Date, x)
%
% [Result, Param] = ini_HydroMod3(Switches, Date, x)
%
% CREC initialization, for details see associated CREC.m file
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

% Fixed parameter
Param.XF = 245;

% Initialization of the reservoir states
%
Param.S = 250;
Param.R = 10;
Param.T = 100;

% Apply CREC for all time steps
%
lP    = length( Date ) ;
Result.Qs    = zeros( lP,1 ) ;

if Switches.exportLight.on == 0
    % Flow components
    Result.inter = zeros( lP, 9 ) ;
    Result.interq = zeros( lP, 15 ) ;
end