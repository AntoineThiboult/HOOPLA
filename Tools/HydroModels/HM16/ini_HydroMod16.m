function [Result, Param] = ini_HydroMod16(Switches, Date, x)
%
% [Result, Param] = ini_HydroMod16(Switches, Date, x)
%
% HydroMod16 initialization, for details see associated HydroMod16.m file
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

% Fixed parameters
Param.N = 16;
Param.XF1 = 25;
Param.XF2 = 200;

% Initialization of the reservoir states
%
Param.S = 100;
Param.L = 50;
Param.T = 20;
Param.s1 = 2;
Param.s2 = 2;
Param.s3 = 2;
Param.s4 = 2;
Param.s5 = 2;
Param.s6 = 2;
Param.s7 = 2;
Param.s8 = 2;
Param.s9 = 2;
Param.s10 = 2;
Param.s11 = 2;
Param.s12 = 2;
Param.s13 = 2;
Param.s14 = 2;
Param.s15 = 2;
Param.s16 = 2;

% Initialize HydroMod16 for all time steps
%
lP    = length( Date ) ;
Result.Qs    = zeros( lP,1 ) ;

if Switches.exportLight.on == 0
    % Flow components
    Result.inter = zeros( lP, 8 ) ;
    Result.interq = zeros( lP, 15 ) ;
end
