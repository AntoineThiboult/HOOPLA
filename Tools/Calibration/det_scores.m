function [perf, varargout]  = det_scores( Q, Qs )
% 
% [perf, varargout]  = det_scores( Q, Qs )
% 
% Compute deterministic score for calibration purpose
% 
% Inputs : 
%   - Q  : observed streamflow (target)
%   - Qs : simulated streamflow
% 
% Output : 
%   - perf : 1x17 matrix containing the following scores
%                    1    2        3       4   5       6      7   8   9       10     11  12     13      14    15    16     17       
%           perf = [RMSE RMSEsqrt RMSElog MSE MSEsqrt MSElog MAE NSE NSEsqrt NSEinv PVE PVEabs Balance r     bKGE  gKGE   KGEm] ;
%   - varargout(1) = modified version of the scores where the optimal result is 0. ex NSE -> abs(NSE-1)
% 
% Programmed by G. Seiller. 
% Revision: 0.0 Date: 2017 Antoine Thiboult
%   - Correction formulation PVEabs
%   - Formating

%% Initialization
ine0 = find(Q >= 0); % ind not empty
ineN = find(isnan(Q)==0); % check for NaNs
ines0 = find(Qs >= 0); % ind not empty
inesN = find(isnan(Qs)==0); % check for NaNs
ine = sort(intersect(intersect(ine0,ineN),intersect(ines0,inesN)));

Qobs = Q(ine,:) ; % Qobs
Qsim = Qs(ine,:) ; % QSim

%% Scores
% Error
RMSE = sqrt(sum((Qsim-Qobs).^2)/length(Qobs)); % Root Mean Square Error
RMSEsqrt = sqrt(sum((sqrt(Qsim)-sqrt(Qobs)).^2)/length(Qobs)); % Root Mean Square Error on sqrt(Q)
RMSElog = sqrt(sum((log(Qsim+eps)-log(Qobs+eps)).^2)/length(Qobs)); % Root Mean Square Error on sqrt(Q)
MSE  = mean((Qobs-Qsim).^2); % Mean Square Error
MSEsqrt = mean((sqrt(Qobs)-sqrt(Qsim)).^2); % Mean Square Error on sqrt(Q)
MSElog = mean((log(Qobs+eps)-log(Qsim+eps)).^2); % Mean Square Error on log(Q)
MAE  = mean(abs(Qobs-Qsim)); % Mean Absolute Error

% Nash-Sutcliffe
NSE = 1-(sum((Qobs-Qsim).^2)/sum((Qobs-mean(Qobs)).^2)); % Nash Sutcliffe Efficiency
NSEsqrt = 1-(sum((sqrt(Qobs)-sqrt(Qsim)).^2)/sum((sqrt(Qobs)-mean(sqrt(Qobs))).^2)); % Nash Sutcliffe Efficiency on sqrt(Q)
NSEinv = 1-(sum(((1./(Qobs+mean(Qobs)/100))-(1./(Qsim+mean(Qobs)/100))).^2)/sum(((1./(Qobs+mean(Qobs)/100))-mean((1./(Qobs+mean(Qobs)/100)))).^2)); % Nash Sutcliffe Efficiency on inverse(Q)

% Water Balance
PVE = (sum(Qsim-Qobs)/sum(Qobs))*100; % Percentage Volume Error (%)
PVEabs = (sum(abs(Qsim-Qobs))/sum(Qobs))*100; % Absolute Percentage Volume Error (%)
Balance = 1-abs(sqrt(sum(Qsim)/sum(Qobs))-sqrt(sum(Qobs)/sum(Qsim))); % Water balance

% Modified Kling-Gupta Efficiency
r = corrcoef(Qobs,Qsim); r = r(2,1); % Linear correlation coefficient (best 1)
bKGE = mean(Qsim)/mean(Qobs); % Beta KGE : Bias (best 1)
gKGE = (std(Qsim)/mean(Qsim))/(std(Qobs)/mean(Qobs)); % Gamma KGE : Variation coefficient ratio (best 1)
KGEm = 1-sqrt(((r-1)^2)+((bKGE-1)^2)+((gKGE-1)^2)); % Modified Kling-Gupta Efficiency (best 1)

%% Performance file
%       1    2        3       4   5       6      7   8   9       10     11  12     13      14    15    16     17       
perf = [RMSE RMSEsqrt RMSElog MSE MSEsqrt MSElog MAE NSE NSEsqrt NSEinv PVE PVEabs Balance r     bKGE  gKGE   KGEm] ;
if nargout == 2 
    varargout(1) = {abs(perf - [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0 ,0, 1, 1, 1, 1])}; % perf not oriented with objective 0
end