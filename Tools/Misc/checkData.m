function [DataObs, varargout] = checkData(Switches, DataPath)
%
% [dataObs, DataFcast] = checkData(Switches, DataPath)
%
% Check data availability. If not provided and essential, throw and error
%
% Inputs :
%   Switches : structure containing information about
%      calibration/simulation/forecasting preferences
%   DataPath : Path to the data
%
% Outputs :
%   DataObs : Structure containing all the observation data necessary for
%      calibration/simulation.
%   DataFcast : Structure containing all the observation data necessary for
%      forecasting.
%
% Programmed by A. Thiboult (2016)

%% Check for observed Data
if Switches.verb.on; dispstat('Verifying that all necessary data are provided...','keepprev');end
runMode=dbstack; % Check if checkData is called from ini_cal/sim/fcast
DataObs=load(DataPath.dataObs);

% General
if isfield(DataObs,'Date')==0;             error('Hydrology:Data',strcat('\nDates not provided\n'));                        end
if isfield(DataObs,'Pt')==0;               error('Hydrology:Data',strcat('\nPrecipitations not provided\n'));               end

% Calibration
if Switches.calibration.on == 1
    if isfield(DataObs,'Q')==0;            error('Q data not provided. No calibration possible - error');                   end
else
    if isfield(DataObs,'Q')==0         
        warning('Hydrology:Data','Q data not provided, set to NaN'); 
        DataObs.Q= NaN(length(DataObs.Date),1); 
    end
end

% Potential evapotranspiration
if Switches.petCompute.on == 1
    listPET = find(Switches.isE)';
    for iPET = listPET
        listVarPET = strsplit( sprintf('%s_%s',Switches.nameE{iPET,2:3}), '_' );  % Necessary data for the iPET'th PET formulation
        for iVarPET = 1 : numel(listVarPET)
            if isfield( DataObs, listVarPET(iVarPET) ) == 0;
                error('PET:Data','\n%s data not provided for %s PET - error\n', listVarPET{iVarPET}, Switches.nameE{iPET,1} );
            end
        end
    end
    if isfield(DataObs,'E')==0;            DataObs.E = NaN(length(DataObs.Date),1);                                         end
elseif Switches.petCompute.on == 0
    if isfield(DataObs,'E')==0;            error('Hydrology:Data',strcat('\nPotential evapotranspiration not provided\n')); end
end

% Snowmelt accounting
if Switches.snowmeltCompute.on == 1
    listSAR = find(Switches.isS)';
    for iSAR = listSAR
        listVarSAR = strsplit( sprintf('%s_%s',Switches.nameS{iSAR,2:3}), '_' );  % Necessary data for the iSAR'th SAR formulation
        for iVarSAR = 1 : numel(listVarSAR)
            if isfield( DataObs, listVarSAR(iVarSAR) ) == 0;
                error('SAR:Data','\n%s data not provided for %s SAR - error\n', listVarSAR{iVarSAR}, Switches.nameS{iSAR,1} );
            end
        end
    end
    
    if isfield(DataObs,'Tmin')==0 && strcmp(Switches.nameS(Switches.isS),'CemaNeige')
        DataObs.Tmin=NaN(size(DataObs.Date,1),1);
        warning('Hydrology:Data','\n%s % s %s\n',...
            'Tmin not provided. Tmin set to NaN.',...
            'CemaNeige: Because Tmin is missing, the USGS function is used to compute snow fraction.');
    end
    if isfield(DataObs,'Tmax')==0 && strcmp(Switches.nameS(Switches.isS),'CemaNeige')
        DataObs.Tmax=NaN(size(DataObs.Date,1),1);
        warning('Hydrology:Data','\n%s % s %s\n',...
            'Tmax not provided. Tmax set to NaN.',...
            'CemaNeige: Because Tmax is missing, the USGS function is used to compute snow fraction.');
    end
end

% Meteorological forecast
if strcmp(runMode(2).name,'ini_forecast') && Switches.forecast.perfectFcast.on == 0
    if Switches.forecast.metEns.on == 1
        metFile=matfile(DataPath.dataMetFcast);
        DataMetFcast=metFile.Met_fcast(1,DataPath.dataMetFcastEnsMb);
    else
        DataMetFcast=load(DataPath.dataMetFcast);
    end
    if isfield(DataMetFcast,'Date')==0;                 error('Hydrology:Data',strcat('\nMeteorological date matrix not provided\n'));                end
    if isfield(DataMetFcast,'leadTime')==0;             error('Hydrology:Data',strcat('\nMeteorological lead times not provided\n'));                 end
    if isfield(DataMetFcast,'Pt')==0;                   error('Hydrology:Data',strcat('\nMeteorological precipitation forecast not provided\n'));     end
    if isfield(DataMetFcast,'T')==0;                    error('Hydrology:Data',strcat('\nMeteorological mean temperature not provided\n'));           end
    if isfield(DataMetFcast,'Tmin')==0
        warning('Hydrology:Data',strcat('\nTmin meteorological forecast not provided. Tmin set to NaN.\n'));
        DataMetFcast.Tmin=DataMetFcast.T*NaN;
        if Switches.snowmeltCompute.on == 1 && strcmp(Switches.nameS(Switches.isS),'CemaNeige')
            warning('Hydrology:Data','\n%s %s\n',...
                'Because the meteorological forecast for Tmin is missing, the USGS function is used to compute snow fraction.',...
                'This may result in a decrease of performance, especially if the the Hydrodel function was used during calibration');
        end
    end
    if isfield(DataMetFcast,'Tmax')==0
        warning('Hydrology:Data',strcat('\nTmax meteorological forecast not provided. Tmax set to NaN.\n'));
        DataMetFcast.Tmax=DataMetFcast.T*NaN;
        if Switches.snowmeltCompute.on == 1 && strcmp(Switches.nameS(Switches.isS),'CemaNeige')
            warning('Hydrology:Data','\n%s %s\n',...
                'Because the meteorological forecast for Tmax is missing, the USGS function is used to compute snow fraction.',...
                'This may result in a decrease of performance, especially if the the Hydrodel function was used during calibration');
        end
    end
    if Switches.forecast.hor>size(DataMetFcast.Pt,2);   error('Hydrology:Data',strcat('\nThe specified forecast horizon is longer than the meteorological forecast horizon\n')); end
    varargout{1}=DataMetFcast;
elseif Switches.forecast.on == 1 && Switches.forecast.perfectFcast.on == 1
    varargout{1}=[];
end
