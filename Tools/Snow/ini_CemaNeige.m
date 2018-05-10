function [SarResult, SarParam] = ini_CemaNeige(Switches, Data, x)
%
% [SarParam, SarResult] = CemaNeigeInit(Switches, data, x)
%
% Initialize Cemaneige and preallocates space
% 
% Coded by G. Seiller (2013)
% Slightly modified by A. Thiboult (2017)


%Cemaneige Parameters
SarParam.CTg = x(end-1);
SarParam.Kf = x(end);

SarParam.G     = zeros( 1,5 ) ;
SarParam.eTg   = zeros( 1,5 ) ;

SarParam.Zz     = Data.Zz5;
SarParam.ZmedBV = Data.Zz5(3);
SarParam.Beta   = Data.Beta;
SarParam.gradT  = Data.gradT;
SarParam.Tf     = 0;
SarParam.QNBV   = Data.QNBV;
SarParam.Vmin   = Data.Vmin;

% Preallocation
Nbarea              = length(SarParam.G);
SarResult.runoffD   = zeros(length(Data.Date),1);

if Switches.exportLight.on == 0
    % Preallocation
    SarResult.Pg         = zeros(length(Data.Date),Nbarea);
    SarResult.Pl         = zeros(length(Data.Date),Nbarea);
    SarResult.G          = zeros(length(Data.Date),Nbarea);
    SarResult.eTg        = zeros(length(Data.Date),Nbarea);
    SarResult.snowMelt      = zeros(length(Data.Date),Nbarea);
end
end