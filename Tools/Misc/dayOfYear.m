function [doy]=dayOfYear(iDate)
%
% [doy]=dayOfYear(iDate) 
%
% Compute day of the year (1st jan = 1, 31 dec = 365)
% 
% Input 
%   iDate : date as a 1x6 matrix (yyyy,mm,dd,hh,mm,ss)
% 
% Output
%   doy : day of the year (scalar)
%
% Programmed by A. Thiboult (2016)

doy = datenum(iDate)-datenum([iDate(:,1),ones(size(iDate,1),2)]);
% doy = datenum(iDate) - datenum(yyyy/01/01/00:00:00)
end
