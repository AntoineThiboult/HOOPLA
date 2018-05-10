function [Result, Param] = ini_HydroMod11(Switches, Date, x)
%
% [Result, Param] = ini_HydroMod11(Switches, Date, x)
%
% MORDOR initialization, for details see associated MORDOR.m file
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

% Parameters
Param.x = x ;

% Initialization of the reservoir states
%
Param.U = x(5)*0.5 ;
Param.L = x(6)*0.5 ;
Param.Z = 50 ;
Param.N = 0.5 ;

% UH2 initialization --> An imaginary part is produced when t > 2x4
%
x4c       = ceil( x(4) ) ;
t         = ( 0:x4c )' ;
t2        = ( x4c:ceil(2*x(4)) )' ;
SH2       = [ 0.5*(t(1:end-1)/x(4)).^2.5; 1-0.5*(2-t2(1:end-1)/x(4)).^2.5; 1 ] ;
SH2       = real( SH2 ) ;
Param.UH2 = diff( SH2 ) ;
Param.H2vsal  = zeros( size(Param.UH2) ) ;
Param.H2rur  = zeros( size(Param.UH2) ) ;
Param.H2vn  = zeros( size(Param.UH2) ) ;

% Apply MORDOR for all time steps
%
lP    = length( Date ) ;
Result.Qs    = zeros( lP,1 ) ;

if Switches.exportLight.on == 0
    % Flow components
    Result.inter = zeros( lP, 11 ) ;
    Result.interq = zeros( lP, 15 ) ;
end