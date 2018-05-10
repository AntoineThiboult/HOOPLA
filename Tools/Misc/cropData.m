function [DataObs, varargout] = cropData(Switches, DataObs, iM, iE, iS, varargin)
%
% [DataObs, DataFcast, DataWu] = CropData(Switches, Data, iM, iE, iS, DataFcast)
%
% Perform data croping according to specified calibration/simulation/forecasting
% dates and warm up.
%
% Inputs :
%   Switches : structure containing information about
%      calibration/simulation/forecasting preferences
%   DataObs : Observation data necessary for
%      calibration/simulation/forecast
%   DataObs : Observation data necessary for forecasting
%
% Outputs :
%   DataObs : Structure containing all the observation data necessary for
%      calibration/simulation.
%   DataFcast : Structure containing all the observation data necessary for
%      forecasting.
%   DataWu : Structure containing all the observation data necessary for
%      hydrological warm up.
%
% Programmed by A. Thiboult (2016)
% Modified by A. Thiboult (2017) for varying time step.

%% Cropable data

if Switches.verb.on; dispstat('Data selection...','keepprev');end

listVarHydro = strsplit( Switches.nameM{iM,2}, '_' ); % Necessary data for the hydrological model
if Switches.petCompute.on==1
    listVarPET = strsplit( Switches.nameE{iE,2}, '_' );  % Necessary data for the PET formulation
else
    listVarPET = {};
end
if Switches.snowmeltCompute.on==1
    listVarSar = strsplit( Switches.nameS{iS,2}, '_');  % Necessary data for the SAR formulation
else
    listVarSar = {};
end
cropableDataObs=unique(['Date', 'Q', listVarHydro, listVarPET, listVarSar]);
cropableDataFcast={'Pt','T','Tmax','Tmin'};

%% Dates

% Retrieve calibration/validation/forecast start and end dates
runMode=dbstack;
switch runMode(2).name
    case 'ini_calibration'
        dateStart=Switches.calStart;
        dateEnd=Switches.calEnd;
    case 'ini_simulation'
        dateStart=Switches.simStart;
        dateEnd=Switches.simEnd;
    case 'ini_forecast'
        dateStart=Switches.fcastStart;
        dateEnd=Switches.fcastEnd;
end

% Numerical time step
ts = str2double( regexprep(Switches.timeStep,'h','') );

% Regular dates
numericDate =       datenum(DataObs.Date);
numericDateStart =  datenum(dateStart);
numericDateEnd =    datenum(dateEnd);
selec = ismember(numericDate, numericDateStart : ts/24 : numericDateEnd); % indices of dates of interest

% Dates warm up
numericDateStartWarmUp = numericDateStart-(ts/3)*365;
numericDateEndWarmUp   = numericDateStart-ts/24;
selecWarmUp = ismember(numericDate, numericDateStartWarmUp : ts/24 : numericDateEndWarmUp); % indices of dates of the warm up

% Check for discrepancies in dates
if numericDateStart<min(numericDate) || numericDateStart>max(numericDate) || numericDateEnd>max(numericDate) || numericDateEnd<min(numericDate)
    error('Hydrology:Dates',strcat('\nThe specified calibration/simulation/forecasting dates are out of the available period\n'))
end


%% Croping data for forecast

if strcmp(runMode(2).name,'ini_forecast')
    
    % To find the validity time of DataMetFcast.variable(iDate,iLeadTime) add
    % DataMetFcast.leadTime(iLeadTime) to DataMetFcast.Date(iDate,:)
    % validityTime = datevec( datenum(DataMetFcast.Date(iDate,:)) + DataMetFcast.leadTime(iLeadTime) )
    
    if Switches.forecast.perfectFcast.on == 0  % Regular meteorological forecast
        % Input arguments
        TmpDataMetFcast=varargin{1};
        
        % Forecast dates
        numericDateFcast = datenum(TmpDataMetFcast.Date);
        numericDateStartFcast =  datenum(dateStart);
        numericDateEndFcast   =  datenum(dateEnd);
        
        selecFcast1 = ismember(numericDateFcast, numericDateStartFcast:ts/24:numericDateEndFcast); % indices of dates between forecast start and forecast end
        selecFcast2 = ismember(numericDateFcast-floor(numericDateFcast), Switches.forecast.issueTime./24); % indices of dates corresponding to hydrological issue time
        selecFcast  = selecFcast1 & selecFcast2; % indices that match both conditions (forecasting period and issue time)
        dateFcast  = TmpDataMetFcast.Date(selecFcast,:); % dates for which a meteorological forecast is issued and included within the specified forecasting period
        
        % Error handling
        if sum(selecFcast)==0 && numel(Switches.forecast.issueTime)==1 % throws an error if no meteorological forecast is available at this time
            error('Hydrology:Forecast','\n%s %d. %s\n','There is no available meteorological forecast issued at',...
                Switches.forecast.issueTime)
        end
        if mod(ts,24)~=0 % check meteorological forecast issuing time for hydrological forecast with time steps smaller than 24 hours
            for iIssueTime=1:numel(Switches.forecast.issueTime) % warning if the user ask for a forecast issue time that is not available
                if ~any(dateFcast(:,4) == Switches.forecast.issueTime(iIssueTime))
                    warning('Hydrology:Forecast','\n%s %d. %s\n','There is no available meteorological forecast issued at',...
                        Switches.forecast.issueTime(iIssueTime),'No hydrological forecast will be issued at this time')
                end
            end
        end
        
        dateRef=datevec(numericDate(selec)); % vector of dates containing all times steps during the forecasting period
        selecRef= ismember(datenum(dateRef),datenum(dateFcast)); % find indices of dateRef that also belong to dateFcast
        
        % Initialization
        for iCrop=1:numel(cropableDataFcast)
            DataMetFcast.(cropableDataFcast{iCrop})=NaN(length(dateRef),Switches.forecast.hor);
        end
        
        % Retrieve data
        for iCrop=1:numel(cropableDataFcast)
            DataMetFcast.(cropableDataFcast{iCrop})(selecRef,:)=TmpDataMetFcast.(cropableDataFcast{iCrop})(selecFcast,1:Switches.forecast.hor);
        end
        DataMetFcast.Date=dateRef;
        DataMetFcast.leadTime=TmpDataMetFcast.leadTime(1:Switches.forecast.hor);
        
    elseif Switches.forecast.perfectFcast.on == 1 % Perfect meteorological forcing
        % Creates "fake" meteorological forecast based on observations.
        
        % Forecast dates
        numericDateStartFcast =  datenum(dateStart);
        numericDateEndFcast   =  datenum(dateEnd);
        selecFcast = ismember(numericDate, numericDateStartFcast:ts/24:numericDateEndFcast); % indices of dates of interest
        idSelecFcast = find(selecFcast==1);
        
        if datenum(dateEnd) + Switches.forecast.hor*ts/24 > datenum(DataObs.Date(end,:))
            error('Hydrology:Dates','\n%s %s %s\n','The specified forecasting dates are out of the available period.',...
                'When using perfect forecast, the last forecasting date plus the forecast horizon should not',...
                'exceed the last available day of observed data. Consider changing the end of the forecast period.')
        end
        
        % Initialization of matrices
        DataMetFcast.Pt   = NaN(numel(idSelecFcast),Switches.forecast.hor);
        DataMetFcast.T    = NaN(numel(idSelecFcast),Switches.forecast.hor);
        DataMetFcast.Tmin = NaN(numel(idSelecFcast),Switches.forecast.hor);
        DataMetFcast.Tmax = NaN(numel(idSelecFcast),Switches.forecast.hor);
        
        % Retrieve "forecast" data from observation data
        for iSelec=1:numel(idSelecFcast)
            DataMetFcast.Pt(iSelec,:)   = DataObs.Pt(idSelecFcast(iSelec)+1:idSelecFcast(iSelec)+Switches.forecast.hor,:);
            DataMetFcast.T(iSelec,:)    = DataObs.T(idSelecFcast(iSelec)+1:idSelecFcast(iSelec)+Switches.forecast.hor,:);
            DataMetFcast.Tmax(iSelec,:) = DataObs.Tmax(idSelecFcast(iSelec)+1:idSelecFcast(iSelec)+Switches.forecast.hor,:);
            DataMetFcast.Tmin(iSelec,:) = DataObs.Tmin(idSelecFcast(iSelec)+1:idSelecFcast(iSelec)+Switches.forecast.hor,:);
        end
        DataMetFcast.Date               = DataObs.Date(selecFcast,:);
        DataMetFcast.leadTime           = (1:Switches.forecast.hor)*ts/24;
        
        if mod(ts,24)~=0 % Set to NaN values that corresponds to a date that isn't a date when a hydrological forecast is issued
            % (only when modeling time step is smaller than 24 hours, otherwise 1 forecast per day is issued)
            idNaN = ~ismember(DataMetFcast.Date(:,4),Switches.forecast.issueTime);
            for iCrop=1:numel(cropableDataFcast)
                DataMetFcast.(cropableDataFcast{iCrop})(idNaN,:)=NaN;
            end
        end
    end
    
    % Output arguments
    varargout{1}=DataMetFcast;
end


%% Croping data for Warm up

if Switches.warmUpCompute.on == 1 % Warm up
    
    if sum(selecWarmUp)==365*8    % check if one year (3h time step) or 8 years (24h time step) is available prior to date -- case YES
        
        for iCrop=1:numel(cropableDataObs)
            DataWu.(cropableDataObs{iCrop})=DataObs.(cropableDataObs{iCrop})(selecWarmUp,:);
        end
        DataWu=mergeStruct(DataObs,DataWu);
        
    elseif sum(selecWarmUp)-sum(selec)<365*8 % check if one year (3h time step) or 8 years (24h time step) is available prior to date -- case NO
        
        warning('Hydrology:Continuity','\n%s %s %i %s %s \n%s','No data available before calibration/simulation/forecasting dates',...
            'or available time period preceding date shorter than', ts/3 ,'year (~3000 time steps).',...
            'Calibration/Simulation/Forecasting period data are used to create a  mean year for warm up.',...
            'Please check for hydrological continuity between the end of warm up and the beginning of forecast period.')
        
        % Creation of sliding years which ends by the day of the year that just precedes the beginning of the simulation
        doy = dayOfYear(DataObs.Date); % day of the year of the entire period where catchment data are available
        doyStart=datenum(dateStart)-datenum(strcat(dateStart(1:4),'/01/01/00:00:00'));% day of the year of the first day of simulation
        doyEndYearWu=mod(doyStart-ts/24,365);   % day of the year of the last timestep of the warm up
        idEndWu=find(doy==doyEndYearWu); % indices of the end of each sliding year
        idStartWu=idEndWu-365*24/ts+1;   % indices of the start of each sliding year
        idStartWu(1)=[]; % Removes the 1rst year first timestep index as there is no corresponding 1rst year start timestep
        idEndWu(1)  =[]; % Removes the 1rst year last timestep index as there is no corresponding first timestep
        
        % Find the average sliding year in terms of precipitation
        yearlyPt=NaN(numel(idStartWu),1);
        for iYear=1:numel(idStartWu)
            yearlyPt(iYear)=sum(DataObs.Pt(idStartWu(iYear):idEndWu(iYear))); % compute the total precipitation of each sliding year
        end
        [~,avgYearPt]=min(abs(yearlyPt-mean(yearlyPt)));
        idWu=idStartWu(avgYearPt):idEndWu(avgYearPt);
        
        % Retrieve the data of the sliding average warm up year
        for iCrop=1:numel(cropableDataObs)
            DataWu.(cropableDataObs{iCrop})=DataObs.(cropableDataObs{iCrop})(idWu,:);
        end
        DataWu=mergeStruct(DataObs,DataWu);
        
    end
    
    % Output arguments
    varargout{2}=DataWu;
end


%% Croping observed data for cal/sim/fcast

for iCrop=1:numel(cropableDataObs)
    DataObs.(cropableDataObs{iCrop})=DataObs.(cropableDataObs{iCrop})(selec,:);
end

%% Checking ratio of streamflow NaN
if sum(isnan(DataObs.Q))/numel(DataObs.Q)>0.75
    warning('Hydrology:StreamflowData','\nOnly %.1f %% of streamflow are not NaN. Consider changing the period\n',...
        100*(1-sum(isnan(DataObs.Q))/numel(DataObs.Q)))
    if sum(~isnan(DataObs.Q)) == 0 && strcmp(runMode(2).name,'ini_calibration')==1
        error('All streamflow are NaN. The calibration is not possible. Please, change calibration dates')
    end
end