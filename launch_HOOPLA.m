function launch_HOOPLA
%
% Hoopla documentation is situated in ./Doc. It provides information about
% Hoopla functionalities, input and output formating, etc.
%
% This function is the main function which should be edited to specify
% user's request. Multiple options are available for for the hydrological
% calibration, simulation, and forecast. These options should be specified
% in this file by modifying the values of the structure "Switches".
%
% We did our best to provide you with a bug-free code. Yet,  it is possible
% that some errors may not be picked up. In the case you find a bug, please
% report it to HOOPLA@fsg.ulaval.ca, and attach the launch_Hoopla function
% as well as the error message.
%
% Programmed by A. Thiboult (2017)

clear
close all
clc

% Make tools available
addpath(genpath('Tools'));

%% Section 1: Switches: 1=yes/on, 0=no/off [EDIT THIS SECTION]

% Calibration / Simulation / Forecast
Switches.calibration.on =0; % Run calibration
Switches.simulation.on  =1; % Run simulation
Switches.forecast.on    =0; % Run forecast

% Dates:format 'yyyy/mm/dd/HH:MM:SS' for 3h time step, 'yyyy/mm/dd' for 24h time step
Switches.calStart   = '1997/01/01/03:00:00'; % Beginning of calibration period
Switches.calEnd     = '2007/01/10/00:00:00'; % End of calibration period
Switches.simStart   = '2010/01/01/03:00:00'; % Beginning of simulation period
Switches.simEnd     = '2015/01/01/00:00:00'; % End of simulation period
Switches.fcastStart = '2015/04/15/03:00:00'; % Beginning of forecast period
Switches.fcastEnd   = '2016/07/01/00:00:00'; % End of forecast period

% General switches
Switches.timeStep           ='3h'; % Computation time step. Choose between '3h' and '24h'
Switches.petCompute.on      =1;     % Compute PET
Switches.snowmeltCompute.on =1;     % Compute snowmelt
Switches.warmUpCompute.on   =1;     % Add warm up before modelling
Switches.verb.on            =1;     % Verbose. Display information about computing
Switches.exportLight.on     =1;     % Export fewer data(/results) to save space
Switches.overWrite.on       =1;     % Overwrite existing files created by HOOPLA
Switches.parallelCompute.on =0;     % Parallel computing

% Calibration switches
Switches.calibration.export.on = 1;    % Export calibrated parameters to ./Data for future Simulation/Forecast once calibration is performed
Switches.calibration.snowCal.on = 0;   % Calibrate snow module (if 0, default values are used)
Switches.calibration.method ='SCE';    % Choose between 'DDS' and 'SCE'
Switches.calibration.rmWinter.on = 1;  % Remove the Quebec "ice months" (dec, jan, fev, mar)
Switches.calibration.score ='NSE';     % Performance criteron (RMSE,MSE,NSE,etc. Enter "help det_score" in terminal to see all available scores)
Switches.calibration.maxiter = 50;    % Maximum number of itération during calibration
Switches.calibration.SCE.ngs = 50;     % Number of Complexes for the SCE optimization

% Forecast switches
Switches.forecast.issueTime        =6;    % Hour of the day for which a forecast is issued (can be several per day ex: [6 12 18 24])
Switches.forecast.perfectFcast.on  =0;    % Use meteorological observations as meteorological forecast
Switches.forecast.hor              =80;  % Horizon of the forecast (in time steps)
Switches.forecast.metEns.on        =0;    % Use meteorological ensemble forecast

% Data Assimilation switches
Switches.DA.on=1;               % Perform data assimilation
if Switches.DA.on == 1
    Switches.DA.tech='EnKF';    % Choose either 'EnKF' (Ensemble Kalman Filter), 'PF' (particle filter), or 'PertOnly' (perturbation of inputs only)
    Switches.DA.Uc_Q=0.1;       % Discharge (standard deviation=10% * Qobs)
    Switches.DA.Uc_Pt=0.5;      % Rainfal (standard deviation=50% * Pt mm, reference:Liu et al 2012,Reichl 2002)
    Switches.DA.Uc_Tpet=2;      % temperature for PET (std deg cel)
    Switches.DA.Uc_Tsno=2;      % temperature for snow melt(std in deg cel)
    Switches.DA.Uc_Tmax=2;      % max temperature for snow melt(std in degree cel)
    Switches.DA.Uc_Tmin=2;      % min temperature for snow melt(std in degree cel)
    Switches.DA.Uc_E=0.1;       % PET (standard deviation=10% * E). Value used only if Switches.petCompute.on=0
    Switches.DA.dt=8;           % delta t between two correction steps
    Switches.DA.N=50;          % Ensemble size
    switch Switches.DA.tech
        case 'EnKF'
        case 'PF'
            Switches.DA.PF.ResampTech='systematic_resampling'; % resampling technique. Either 'multinomial_resampling' or 'systematic_resampling'
            Switches.DA.PF.resampleThresh = inf;               % Particle effective ensemble size before resampling. Included in [0 N]. 0 = never resample, ..., N = resample at each time step
        case 'PertOnly'
    end
end

%% Section 2: Loading framework tools [DO NOT EDIT THIS SECTION]
Switches=list_models(Switches);
Switches=list_catchments(Switches);

%% Section 3: Catchment, hydrological model, PET model, and snow model selection [EDIT THIS SECTION]

% Catchments
Switches.isC=zeros(size(Switches.nameC,1),1); % Switch for the catchments. If isC(n)==1, catchment n is ran
Switches.isC(1)=1;

% Hydrological models
Switches.isM=zeros(size(Switches.nameM,1),1); % Switch for the models. If isM(n)==1, model n is ran
Switches.isM(:)=1;

% PET models
Switches.isE=zeros(size(Switches.nameE,1),1); % Switch for the models. If isE(n)==1, model n is ran
Switches.isE(1)=1;

% Snow models
Switches.isS=zeros(size(Switches.nameS,1),1); % Switch for the models. If isS(n)==1, model n is ran
Switches.isS(1)=1;

%% Section 4: Run HOOPLA [DO NOT EDIT THIS SECTION]
run_HOOPLA(Switches)

