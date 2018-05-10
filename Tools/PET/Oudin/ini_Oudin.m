function [PetData] = ini_Oudin(Switches, Data)
% [PetData] = ini_Oudin(Data)
%
% Compute the prerequisite values that are necessary for the computation of Oudin
%
%   Inputs
%       Data.
%           Date   = Date / date (datenum)
%           T      = Temperature moyenne journaliere / daily mean temperature (°C)
%           Lat    = Latitude de la station / station latitude (°)
%
%   Outputs
%       PetData.
%           Re     = Extraterrestrial radiation (MJ/m2/time unit)
%           lambda = Latent vaporization energy (MJ/kg)
%           rho    = Water volumic mass (= 1000 kg/L)
%           T      = Mean air temperature
%           DL     = Maximum day light (h)
%
% Coded by A. Thiboult (2017)

%% Computation day of year (1-365)
doy = dayOfYear(Data.Date);

%% Constants
Gsc = 0.082; % (MJ/m2/min)
rho = 1000 ; % (kg/L)
lambda = 2.501-0.002361.*Data.T; % (MJ/kg)

%% Computation
Lrad = (pi()*Data.Lat)/180; % (rad)
ds = 0.409*sin((2*pi()/365)*doy-1.39); % (rad)
dr = 1+0.033*cos(doy*2*pi()/365); % no unit
omega = acos(-tan(Lrad).*tan(ds)); % (rad)

switch Switches.timeStep
    case '3h'
        b = 2.*pi().* (doy-81)./ 364; % component of the seasonal correction.
        Sc = 0.1645 * sin(2* b) - 0.1255 * cos(b) - 0.025 * sin(b); % seasonal correction
        t = Data.Date(:,4)+0.5; % standard clock time at the midpoint of the period [hour]
        Lz = 75; % longitude of the centre of the local time zone [degrees west of Greenwich]
        Lm = 72.0; % longitude of the measurement site [degrees west of Greenwich]
        omega0 = pi/12 .* (t+ 0.06667.*(Lz-Lm)+ Sc - 12); % solar time angle at midpoint of the period
        t1 = 3; % time step (in hours)
        omega1 = omega0 - (pi * t1)/24; % solar time angles at the beginning of the period
        omega2 = omega0 + (pi * t1)/24; % solar time angles at the end of the period
        Re = (12*60/pi).*Gsc.*dr.*...
            ((omega2-omega1).*sin(Lrad).*sin(ds)+cos(Lrad).*cos(ds).*(sin(omega2)-sin(omega1))); % (3-hour) solar radiation formula from Allen et al.(1998) corrected by A. Thiboult
        
    case '24h'
        Re = (24*60/pi())*Gsc.*dr.*(omega*sin(Lrad).*sin(ds)+cos(Lrad).*cos(ds).*sin(omega)); % in MJ/m2/day
end


%% Outputs
PetData.Re     = Re;
PetData.lambda = lambda;
PetData.rho    = rho;
PetData.T      = Data.T;

end

