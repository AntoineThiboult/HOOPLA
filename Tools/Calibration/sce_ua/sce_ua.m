function [bestx,bestf,allbestf] = sce_ua(CostFunctionName,~,VarMin,VarMax,nComplex,userdata)
%        [bestx,bestf,allbest] = sce_ua(CostFunctionName,~,VarMin,VarMax,nComplex,userdata)
%
% Shuffled Complex Evolution (SCE-UA) METHOD
%
%
% INPUTS
%  fctname = character string of the function to optimize
%  x0      = the initial parameter array at the start;
%          = the optimized parameter array at the end;
%  VarMin  = the lower bound of the parameters
%  VarMax  = the upper bound of the parameters
%  nComplex = number of complexes (sub-pop.)- between 2 and 20
%  userdata (optional) - for the function to optimize
%
% OUTPUTS
%  bestx   = best parameters
%  bestf   = best value of the cost function
%  allbest = best value of the cost function for each iteration
%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license_sce.txt" for license terms.
% Project Code: YPEA110
% Project Title: Implementation of Shuffled Complex Evolution (SCE-UA)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%
% Modification: 20/09/2017 by Antoine Thiboult.
%               - modification of the file sce_ua for compatibility with ls
%                 HOOPLA
%               - add a convergency criteria condition

%% Problem Definition

nVar = size(VarMax,2);      % Number of Unknown Variables
VarSize = [1 nVar];         % Unknown Variables Matrix Size
CostFunction = str2func( CostFunctionName );
kStep = 10; % number of steps considered to compute gain for convergency
epsCost = 1e-4; % Minimum improvement of the cost function in the last kStep

%% SCE-UA Parameters

MaxIt = userdata.Switches.calibration.maxiter;  % Maximum Number of Iterations
nPopComplex = 5;                                % Complex Size
nPopComplex = max(nPopComplex, nVar+1);         % Nelder-Mead Standard
nPop = nComplex*nPopComplex;                    % Population Size

I = reshape(1:nPop, nComplex, []);

% CCE Parameters
cce_params.q = max(round(0.5*nPopComplex),2);   % Number of Parents
cce_params.alpha = 3;   % Number of Offsprings
cce_params.beta = 5;    % Maximum Number of Iterations
cce_params.CostFunction = CostFunction;
cce_params.VarMin = VarMin;
cce_params.VarMax = VarMax;
cce_params.userdata = userdata; 

%% Initialization

% Empty Individual Template
empty_individual.Position = [];
empty_individual.Cost = [];

% Initialize Population Array
pop = repmat(empty_individual, nPop, 1);

% Initialize Population Members
for i=1:nPop
    pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
    pop(i).Cost = CostFunction(pop(i).Position, userdata);
end

% Sort Population
pop = SortPopulation(pop);

% Update Best Solution Ever Found
BestSol = pop(1);

% Initialize Best Costs Record Array
BestCosts = nan(MaxIt, 1);

%% SCE-UA Main Loop

for it = 1:MaxIt
    
    % Initialize Complexes Array
    Complex = cell(nComplex, 1);
    
    % Form Complexes and Run CCE
    for j = 1:nComplex
        % Complex Formation
        Complex{j} = pop(I(j,:));
        
        % Run CCE
        Complex{j} = RunCCE(Complex{j}, cce_params);
        
        % Insert Updated Complex into Population
        pop(I(j,:)) = Complex{j};
    end
    
    % Sort Population
    pop = SortPopulation(pop);
    
    % Update Best Solution Ever Found
    BestSol = pop(1);
    
    % Store Best Cost Ever Found
    BestCosts(it) = BestSol.Cost;
    
    % Show Iteration Information
    fprintf('Iteration %i : Best Cost = %0.4f\n',...
        it, abs(userdata.orientScore(userdata.idScore)-BestCosts(it)));
    
    % Check for cost convergency
    if it > kStep
        gainBestCosts = (BestCosts(it) - BestCosts(it-kStep)) ./  BestCosts(it-kStep); 
        if abs(gainBestCosts) < epsCost
            break
        end
    end
    
end

%% Results
bestx = BestSol.Position; 
bestf = abs(userdata.orientScore(userdata.idScore) - BestSol.Cost); 
allbestf = abs(userdata.orientScore(userdata.idScore) - BestCosts); 
