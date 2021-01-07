function [param]=EnKF_func(param, Qs, ~, Qrp, eQ, DA, ~)
%
% [param]=EnKF_func(param, Qs, ~, Q, eQ, DA, ~)
%
% Assimilation of streamflow to update states variables with the Ensemble
% Kalman Filter
%
% Inputs:
%   param : state variables
%   Qs  : simulated streamflow (from prediction step)
%   Qrp : streamflow with random perturbations
%   eQ  : streamflow Gaussian noise
%   DA  : structure containing EnKF specs
%
% Outputs:
%   param: updated state variables
%
% Programmed by A. Thiboult and Philipp Meier (2014) following Mandel, J.,
% Efficient Implementation of the Ensemble Kalman Filter, Report,
% Univ of Colorado, 2006.

%% Ensemble Kalman Filter

% Number members
N=DA.N;

% Initialization of states matrix
X=NaN(numel(DA.updatedRes),DA.N);

% Assignment of states
for i=1:numel(DA.updatedRes)
    X(i,:)=cat(2,param.(DA.updatedRes{i}));
end

% Observation matrix
z=Qs*ones(N,1);

% EnKF computation
HA=Qs-1/N*z*ones(1,N);
Y=Qrp-Qs;
L=eQ*eQ'/(N-1);
P=L+1/(N-1)*HA*(HA)';
M=P^-1*Y;
Z=HA'*M;
A=X-1/N*(X*ones(N,1))*ones(1,N);
Xa=X+1/(N-1)*A*Z;

% Correction of nonsense state values
Xa(Xa<0)=0;

% Set new states variable
for k=1:DA.N
    for i=1:numel(DA.updatedRes)
        param(k).(DA.updatedRes{i})=Xa(i,k);
    end
end

end
