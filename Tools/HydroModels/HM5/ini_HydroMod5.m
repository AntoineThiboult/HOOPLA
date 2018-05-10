function [Result, Param] = ini_HydroMod5(Switches, Date, x)
%
% [Result, Param] = ini_HydroMod5(Switches, Date, x)
%
% GR4H initialization, for details see associated GR4J.m file
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

% Parameters
Param.x = x' ;

% Initialization of the reservoir states
%
Param.B = 0.9 ;
Param.S = 0.5 * Param.x(1) ;
Param.R = 0.5 * Param.x(3) ;

% UH1 initialization --> the last one always has a value of 1
%
x4c       = ceil( x(4) ) ;
t         = ( 0:x4c )' ;
SH1       = [( t(1:end-1)/x(4) ) .^1.25; 1] ;
Param.UH1 = diff( SH1 ) ;
Param.H1  = zeros( size(Param.UH1) ) ;

% UH2 initialization --> An imaginary part is produced when t > 2x4
%
t2        = ( x4c:ceil(2*x(4)) )' ;
SH2       = [ 0.5*(t(1:end-1)/x(4)).^1.25; 1-0.5*(2-t2(1:end-1)/x(4)).^1.25; 1 ] ;
SH2       = real( SH2 ) ;
Param.UH2 = diff( SH2 ) ;
Param.H2  = zeros( size(Param.UH2) ) ;

% Apply GR4J for all time steps
%
lP    = length( Date ) ;
Result.Qs    = zeros( lP,1 ) ;

if Switches.exportLight.on == 0
    % Flow components
    Result.inter = zeros( lP, 9 ) ;
    Result.interq = zeros( lP, 15 ) ;
end