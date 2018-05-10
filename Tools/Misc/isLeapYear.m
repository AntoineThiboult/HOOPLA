function [ily]=isLeapYear(iYear)
%
% [ily]=isLeapYear(iYear)
% 
% Check if the years are leap years
%
% Input : 
%   - iYear : matrix (1xN)
% Output : 
%   - ily : matrix of logical (1xN). 
%       True = leap year, False = not leap year
% 
% Programmed by A. Thiboult (2016)

ily = mod(iYear,4)==0 & mod(iYear,100)~=0 | mod(iYear,400)==0;
 
end