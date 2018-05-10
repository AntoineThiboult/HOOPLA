function [Param, w]=PF_func(Param, Qs, Q, ~, ~, DA, w)
%
% [Param, w]=PF_func(Param, Qs, Q, ~, ~, DA, w)
%
% Assimilation of streamflow to update states variables with the Particle
% Filter (Sequential Importance Resampling)
%
% Inputs
%   param : state variables
%   Qs  : simulated streamflow (from prediction step)
%   Q   : observed streamflow
%   DA  : structure containing particle filter specs
%   w   : particle weights
%
% Outputs
%   Param : updated state variables
%   w : particle weights
%
% Programmed by A. Thiboult based on :
%   - Arulampalam et. al. (2002). A tutorial on particle filters for online
% nonlinear/non-gaussian bayesian tracking. IEEE Transactions on signal
% processing, Vol. 50, No. 2
%   - A Particle Filter tutorial made by Diego Andrés Alvarez Marín
% available at: https://www.mathworks.com/matlabcentral/fileexchange/
% 35468-particle-filter-tutorial?focused=5242709&tab=function

%% Particle filter

% Initialization of matrices
x = NaN(DA.N, numel(DA.updatedRes));
% x : sample for q(x{t}|x{t-1},z{t}). => importance density =~ prior (Eq. 62)
for i=1:numel(DA.updatedRes)
    x(:,i) = cat(1,Param.(DA.updatedRes{i}));
end
% w : weights w{t} <=> w{t-1} * p(z{t}|x{t}) (Eq. 63)
if sum(normpdf(Q-Qs, 0, Q*DA.Uc_Q))~=0
    w = w .* normpdf(Q-Qs, 0, Q*DA.Uc_Q);
else % Allow the simulation to continue in an open loop fashion if all probabilities are null
    w = repmat(1/DA.N, 1, DA.N); % All particles have the same weight
    return
end

% Normalize weights
w = w./sum(w);

%% Resampling

% Effective sample size (Eq. 50)
Neff = 1/sum(w.^2);
% Resampling threshold
Nthresh = DA.PF.resampleThresh*DA.N;
% Resampling
if Neff < Nthresh
    [x, w] = resampleParticle(x, w, DA.PF.ResampTech);
end

%% Set new states variable
for i=1:numel(DA.updatedRes)
    for k=1:DA.N
        Param(k).(DA.updatedRes{i})=x(k, i);
    end
end
end

%% Resampling function
function [x, w, id] = resampleParticle(x, w, resampTech)
%
% [x, w, id] = resampleParticle(x, w, resampTech)
%
% Programmed by A. Thiboult heavily based on 
%   - A Particle Filter tutorial made by Diego Andrés Alvarez Marín
%       available at: https://www.mathworks.com/matlabcentral/fileexchange/
%       35468-particle-filter-tutorial?focused=5242709&tab=function
% and on 
%   - Comparison of Resampling Schemes for Particle Filtering (2005) 
%       Image and Signal Processing and Analysis

N = length(w); % Number of particles

switch resampTech
    case 'multinomial_resampling'
        id = randsample(1:N, N, true, w);
    case 'systematic_resampling'
        histEdges = min([0 cumsum(w)],1); % protect against accumulated round-off
        histEdges(end) = 1; % get the upper edge exact
        u1 = rand/N;
        [~, id] = histc(u1:1/N:1, histEdges);
    otherwise
        error('Unknownd resampling strategy')
end;

x = x(id,:); % extract new particles
w = repmat(1/N, 1, N); % All particles have the same weight

end
