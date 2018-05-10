function [Data, w]=ini_DA(Switches,Data)
%
% [Data, w]=ini_DA(Switches,Data)
%
% The function produces noisy imputs (random perturbations RP) for the EnKF
%
% Inputs:
%   Switches : structure containing fields
%     Switches.DA.on          Switch on/off
%     Switches.DA.Uc_Q        Discharge standard deviation
%     Switches.DA.Uc_Pt       Rainfal standard deviation
%     Switches.DA.Uc_Tpet     Temperature for PET (std deg cel)
%     Switches.DA.Uc_Tsno     Temperature for snow melt(std in deg cel)
%     Switches.DA.Uc_Tmax     Max temperature for snow melt(std in degre)
%     Switches.DA.Uc_Tmin     Min temperature for snow melt(std in degre)
%     Switches.DA.N           Ensemble size
%
% Outputs:
%   Data:
%     Data.Pt   : Observed rainfall
%     Data.T    : Observed mean temperature
%     Data.Tmax : Observed max temperature
%     Data.Tmin : Observed min temperature
%     Data.Q    : Observed streamflow
%     w : initial particles weights
%
% Programmed by A. Thiboult (2016)

%% Uncertainties of data for Data assimilation

% Temperatures
Data.TpetRP = normrnd(repmat(Data.T,1,Switches.DA.N),Switches.DA.Uc_Tpet);
Data.TsnoRP = normrnd(repmat(Data.T,1,Switches.DA.N),Switches.DA.Uc_Tsno);
Data.TmaxRP = normrnd(repmat(Data.Tmax,1,Switches.DA.N),Switches.DA.Uc_Tmax);
Data.TminRP = normrnd(repmat(Data.Tmin,1,Switches.DA.N),Switches.DA.Uc_Tmin);

% Streamflow
Data.QRP    = normrnd(repmat(Data.Q,1,Switches.DA.N),Switches.DA.Uc_Q.*repmat(Data.Q,1,Switches.DA.N)); % Negative strreamflow are possible.
Data.eQRP   = repmat(Data.Q,1,Switches.DA.N)-Data.QRP;

% Rainfall
k = (Data.Pt).^2 ./ (Switches.DA.Uc_Pt.*Data.Pt).^2; % Scale parameter
teta = (Switches.DA.Uc_Pt*Data.Pt).^2 ./ Data.Pt; % Shape parameter

Data.PtRP = gamrnd(repmat(k,1,Switches.DA.N),...
    repmat(teta,1,Switches.DA.N));
Data.PtRP(isnan(Data.PtRP))=0;

% Initialize weights for particle filter
w = ones(1,Switches.DA.N) ./ Switches.DA.N;

% Potential evapotranspiration
if Switches.petCompute.on == 0
    Data.ERP = normrnd(repmat(Data.E,1,Switches.DA.N),Switches.DA.Uc_E);
end

end
