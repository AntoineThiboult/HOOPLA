function [Result, Param] = ini_HydroMod10(Switches, Date, x)
%
% [Result, Param] = ini_HydroMod10(Switches, Date, x)
%
% MOHYSE initialization, for details see associated MOHYSE.m file
%
% Programmed by G. Seiller, Univ. Laval (05-2013)
% Slightly modified by A. Thiboult (2016)

% Parameters
Param.x = x ;
k = 80; % memory

% Initialization of the reservoir states
%
Param.S = 40 ;
Param.R = 30 ;

% UH initialization
%
t        = ( 1:k )' ;
SH       = t.^(x(6)-1).*exp(-t/x(7));
sumSH    = sum(SH);
Param.UH = SH./sumSH;
Param.HUQ1 = zeros( size(Param.UH) ) ;
Param.HUQ2 = zeros( size(Param.UH) ) ;
Param.HUQ3 = zeros( size(Param.UH) ) ;

% Apply MOHYSE for all time steps
%
lP    = length( Date ) ;
Result.Qs    = zeros( lP,1 ) ;

if Switches.exportLight.on == 0
    % Flow components
    Result.inter = zeros( lP, 8 ) ;
    Result.interq = zeros( lP, 15 ) ;
end