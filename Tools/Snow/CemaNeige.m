function [runoffD, SarParam, varargout] = CemaNeige( Pt, T, Tmax, Tmin, Date, SarParam )
%
% [runoffD, SarParam, Pg, Pl, G, eTg, snowMelt] = CemaNeige( Pt, T, Tmax, Tmin, Date, SarParam )
%
% Snow accounting routine. Compute accumulation and snow melt.
%
% Inputs:
%   Pt   = total precipitation (solid + liquid)
%   T    = mean temperature (°C)
%   Tmax = max temperature (°C)
%   Tmin = min temperature (°C)
%   Date = yyyy / mm / dd / hh : mm: ss (1x6 matrix)
%   SarParam.
%       CTg = snow cover thermal coefficient (calibrated paramter)
%       Kf  = snowmelt factor (mm / °C * pdt) (calibrated paramter)
%       G   = snow stock
%       eTg = snow cover thermal state
%       Zz  = catchment elevation quantiles (m)
%       ZmedBV  = median catchment elevation (m)
%       Beta    = correction of precipitation according to elevation (m-1)
%       gradT   = temperature gradient (°C/100 m)
%       Tf      = melting temperature (°C)
%       QNBV    = average annual snow accumulation (mm)
%       Vmin    = percentage of Kf that corresponds to the minimal melting
%                 rate (=0.1 if G = 0)
%
% Outputs: 
%   runoffD = depth of runoff
%   CemaParm = see input above
%   Pg
%   Pl = liquid precipitation
%   snowMelt = snow melt
%
%
% Default parameter values for daily time step (from A. Valery thesis)
%   - CTg : values ranging from 0 to 1 , median = 0.25
%   - Kf : values ranging from 0 to 20 , median = 3.74
% Default parameter values for 3h time step (from Hoopla testing)
%   - CTg : values ranging from 0.75 to 1 , median = 0.93
%   - Kf : values ranging from 0 to 1 , median = 0.40
% Misc: 
%   Gthreshold : quantity of snow above which all the catchment surface is supposed to be covered (mm) 
%            Set to 9/10th of average annual snow accumulation = QNBV*0.9
%
% Following CEMAGREF - A. Valéry (2010)
% 
% Coded by G. Seiller (2013)
% Slightly modified by A. Thiboult (2017)
% Modified by S.Lachance-Cloutier (2022)
%% CemaNeige running

% Variables
G       = SarParam.G;
eTg     = SarParam.eTg;
Zz      = SarParam.Zz;
ZmedBV  = SarParam.ZmedBV;
Beta    = SarParam.Beta;
gradT   = SarParam.gradT;
Tf      = SarParam.Tf;
QNBV    = SarParam.QNBV;
Vmin    = SarParam.Vmin;

% Free parameters
CTg = SarParam.CTg;
Kf = SarParam.Kf;

% Number of elevation bands
nbzalt = 5;

% If it is a leap year, julian days after the 29/02 are shifted by one day for gradT
jd = floor(dayOfYear(Date)); % day of year
if isLeapYear(Date(1))==1 && jd>59
    jd=jd-1;
end
theta = gradT(jd+1);

% Effective temperature
Tz = T + theta*(Zz - ZmedBV)./100;
Tzmax = Tmax+theta*(Zz-ZmedBV)./100;
Tzmin = Tmin+theta*(Zz-ZmedBV)./100;

% Distribution of precipitation over the nbzalt bands
modc = exp(Beta*(Zz-ZmedBV));
c = sum(modc)/nbzalt;
Pz = (1/c)*Pt(1)*exp(Beta*(Zz-ZmedBV));

% Snow franction
Fracneige = zeros(1,nbzalt);
for z = 1 : nbzalt
    if Zz(z) < 1500 && ~any(isnan([Tzmax,Tzmin]))% Hydrotel function
        if Tzmax(z) <= 0
            Fracneige(z) = 1;
        elseif Tzmin(z) >= 0
            Fracneige(z) = 0;
        else
            Fracneige(z) = 1-(Tzmax(z)/(Tzmax(z)-Tzmin(z)));
        end
    else % USGS function (USGS is chosen if Tmax and Tmin are not defined)
        if Tz(z) > 3
            Fracneige(z) = 0;
        elseif Tz(z) < -1
            Fracneige(z) = 1;
        else
            Fracneige(z) = 1-((Tz(z)-(-1))/(3-(-1)));
        end
    end
end
Fracneige = min(Fracneige,1);
Fracneige = max(Fracneige,0);

% Dispatching according to precipitation type
Pg = Pz .* Fracneige;
Pl = Pz - Pg;

% Snow pack updating
G = G + Pg;

% Snow pack thermal state
eTg = CTg*eTg+(1-CTg)*Tz;
eTg = min(0,eTg);

% Melting factor according to snowpack thermal state
fTg = (eTg >= Tf);

% Potential superficial melting
Fpot = (Tz > 0) .* (min(G,Kf*(Tz-Tf).*fTg));

% Ratio of the area covered by snow
Gthreshold = QNBV*0.9;
fnts = min(G/Gthreshold,1);

% Effective melting
snowMelt = Fpot .* ((1-Vmin)*fnts+Vmin);

% Updating of snow stock
G = G - snowMelt;

% Depth total of runoff (sent to the hydrological model)
runoffD = (sum(Pl) + sum(snowMelt))/nbzalt;

% Snowpack state saving
SarParam.G = G;
SarParam.eTg = eTg;

if nargout>2
    % Ouputs
    varargout{1}=Pg;
    varargout{2}=Pl;
    varargout{3}=G;
    varargout{4}=eTg;
    varargout{5}=snowMelt;
end
end