function [varargout] = dds(objfunc_name, sinitial, S_min, S_max, maxiter, userdata)
% function [best, allbest, solution] = dds(objfunc_name, sinitial, S_min, S_max, maxiter, userdata)
% 
% =================================================================================================
% Dynamically Dimensioned Search algorithm by Bryan Tolson.  Version 1.0.2  Oct 2005
%  Small Jan 2014 edit to call get_objfunc.m file when evaluating objective
%  function.
% 
% Coded as a MINIMZER but variable to_min converts a maximization problem
% to a minimization problem without modifying the code below.  In other words,
% specify to_max variable correctly and the code here does not need modification
% for MAX or MIN objectives.
% this algorithm records all sampled solutions but only searches based on
% the best solution
% 
% Inputs:
%       - objfunc_name: string with the objective function name without .m extension
%       - sinitial: initial points to evaluate. Empty vector to force a random initial evaluation
%       - S_min: vector of minimum values of parameters to optimise
%       - S_max: vector of maximum values of parameters to optimise
%       - maxiter: maximum number of iterations
%       - userdata: structure with DDS parameters comfiguration and
%       variables to pass to "objfunc_name" if necessary
% 
% Outputs
%       - best: decision variables found as the best set and its function value
%       - allbest: all best solutions over all iterations
%       - solution: matrix with the following structure
%           col 1: iteration number i
%           col 2: best current solution at iteration i
%           col 3: tested solution at iteration i
%           col4 to col(3+parameters to optimise): parameters tested
% =================================================================================================

% 
% TODO: 
%       - Parallel version
%       - Verify maximization problem: in accumulation results
%       (solution matrix) appers to be used twice the to_max factor.
% Slightly modified by A. Thiboult (2017)

if nargin < 6, userdata = []; end
% parallel evaluation of trials
% if isfield(userdata,'parRuns'), parallel = userdata.parRuns; else parallel = 0; end
% its is the number of function evaluations to initialize the DDS algorithm solution
if isfield(userdata,'evalsInit'), its = userdata.evalsInit; else its = max(5,round(0.005*maxiter)); end
% adjust neighborhood size by modifying the perturbation magnitude
if isfield(userdata,'perturbMagn'), fraction1 = userdata.perturbMagn; else fraction1 = 0.20; end
% variable type flag (0 real, 1 discrete)
if isfield(userdata,'varType'), Discrete_flag = userdata.varType; else Discrete_flag = zeros(length(S_min),1); end
% function orientation (1 to minimisation, -1 to maximisation)
if isfield(userdata,'orientation'), to_max = userdata.orientation; else to_max = 1; end
% verbose
if isfield(userdata,'verbose'), verbose = userdata.verbose; else verbose = 0; end
% evolution percentage to print results
if isfield(userdata,'evolPrint'), evolPrint = userdata.evolPrint; else evolPrint = 10; end

[~, num_dec] = size(S_min); % num_dec is the number of decision variables
solution = zeros(maxiter,3+num_dec); % column 1 =iter#, 2 =best f(x), 3 = test f(x),
% column 4 = dec var 1 value for iter#, 5 = dec var 2 value for iter#, 6 =etc...
S_range=S_max-S_min;

if evolPrint < 1; evolPrint = 100*evolPrint; end
evolProg = round(maxiter*evolPrint/100):round(maxiter*evolPrint/100):maxiter;

% force to use initial supplied solution
if ~isempty(sinitial), its = 1; end

% =================================================================================================
% INITIAL SOLUTION
% =================================================================================================
if its>1   % its is the number of function evaluations to initialize the DDS algorithm solution
    % fprintf('     Finding best starting point for trial %g using %g random samples.\n',trial_num,its);
    ileft=maxiter-its; % use this to reduce number of fevals in DDS loop
    if ileft<=0
        error('#Initialization samples >= Max # function evaluations.')
    end
    for i=1:its
        if Discrete_flag == 0   % continuous variable
            stest = S_min+S_range.*rand(1,num_dec); % uniform random samples
        else                    % discrete case
            for j = 1:num_dec
                stest(1,j) = randi([S_min(1,j), S_max(1,j)],1,1);
            end
        end
        %Jtest = to_max*get_objfunc(stest); % get obj function value
        Jtest = to_max*feval(objfunc_name,stest,userdata); % get obj function value
        if i==1
            Jbest=Jtest;
        end
        if Jtest<=Jbest
            Jbest = Jtest;
            sbest = stest;
        end
        solution(i,1)=i;
        solution(i,2)=to_max*Jbest;
        solution(i,3)=to_max*Jtest;
        solution(i,4:3+num_dec)=stest;
    end
else  % know its=1, using a user supplied initial solution.  Calculate obj func value.
    ileft=maxiter-1; % use this to reduce number of fevals in DDS loop
    stest=sinitial;  % get from the inputs
    Jtest = feval(objfunc_name,stest,userdata); % get obj function value
    Jbest=Jtest;
    sbest=stest;
    solution(1,1)=1;
    solution(1,2)=to_max*Jbest;
    solution(1,3)=to_max*Jtest;
    solution(1,4:3+num_dec)=stest;
end

%it_sbest=its; % needed to initialize variable and avoid code failure when small # iterations
% trial_initial=sbest;  % extra variable here to simplify code for tracking initial DDS solution
% =================================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%  Done initialization procedure in DDS.m

for i=1:ileft  % remaining F evals after initialization
    
    %Determine variable selected as neighbour
    Pn=1.0-log(i)/log(ileft); % 1.0-i/ileft; % probability of being selected as neighbour
    dvn_count=0; % counter for how many decision variables vary in neighbour
    stest = sbest;  % define stest initially as current (sbest for greedy)
    randnums=rand(num_dec,1);
    for j=1:num_dec
        if randnums(j)<Pn % then j th DV selected to vary in neighbour
            dvn_count=dvn_count+1;
            new_value = neigh_value_mixed(sbest(j), S_min(j), S_max(j), fraction1, j, Discrete_flag);
            %new_value=neigh_value(sbest(j),S_min(j),S_max(j),fraction1);
            stest(j)=new_value; % change relevant dec var value in stest
        end
    end
    
    if dvn_count==0 % no DVs selected at random, so select ONE
        dec_var=ceil(num_dec*rand(1,1)); % which dec var to modify for neighbour
        new_value = neigh_value_mixed(sbest(dec_var), S_min(dec_var), S_max(dec_var), fraction1, dec_var, Discrete_flag);
        %new_value=neigh_value(sbest(dec_var),S_min(dec_var),S_max(dec_var),fraction1);
        stest(dec_var)=new_value; % change relevant dec var value in stest
    end
    
    % get ojective function value
    % Jtest = to_max*get_objfunc(stest);
    Jtest = to_max*feval(objfunc_name, stest, userdata);
    
    if Jtest<=Jbest
        Jbest = Jtest;
        sbest = stest;
        %it_sbest=i+its; % iteration number best solution found
        
        %%% write new status file so that best sol'n not lost with long
        %%% runs (i.e. SWAT or other models called).  June 05 - BT
        % Comment this part of code out for fast problems!!
        %filenam='status.out';
        %fid = fopen(filenam,'w'); % opens file and discards current contents
        %zzz=to_max*Jbest;
        %fprintf(fid,'Current best objective function value of %12.5f found at iteration %6.0f\n',zzz,i+its);
        %fprintf(fid,'under parameter set below: \n');
        %fprintf(fid,' %e ',sbest);
        %if verbose
        %    fprintf('Current best objective function value of %12.5f found at iteration %6.0f\n',zzz,i+its);
        %end
        %fprintf('under parameter set below: \n');
        %fprintf(' %e ',sbest);
        %fclose(fid);
        %%%
    end
    
    % accumulate results
    solution(i+its,1)=i+its;
    solution(i+its,2)=to_max*Jbest;
    solution(i+its,3)=to_max*Jtest;
    solution(i+its,4:3+num_dec)=stest;
    
    %%%%%%%% WRITE SCREEN OUTPUT - *** user uncomment if desired ***:
    if any(evolProg==(i+its))
        if verbose
            %fprintf('\n\n');
            prog = round(100*(i+its)/maxiter);  %sum((evolProg==(i+its)).*(evolPrint:evolPrint:100));
            fprintf('Evolution Progress > %3d%% iteration: %3d \t best fitness: %4.3f \t iter fitness: %4.3f \n', ...
                prog, i+its, userdata.orientScore(userdata.idScore)-to_max*Jbest, userdata.orientScore(userdata.idScore)-to_max*Jtest); % AT
            %fprintf('Jtest for iteration %3d is %g\n',i+its,to_max*Jtest);
            %fprintf('Jbest for iteration %3d is %g\n',i+its,to_max*Jbest);
            %fprintf('\n\n');
        end
    end
    
end % DDS function loop

allbest = solution(:,2);
% best = [sbest to_max*Jbest]; AT

if verbose
    fprintf('Best solution found has obj function value of %10.4f \n\n',userdata.orientScore(userdata.idScore)-to_max*Jbest);
end

nout = max(nargout,1);

if nout > 0, varargout(1) = {sbest}; end % AT
if nout > 1, varargout(2) = {userdata.orientScore(userdata.idScore)-to_max*Jbest}; end %AT
if nout > 2, varargout(3) = {userdata.orientScore(userdata.idScore)-allbest}; end
if nout > 3, varargout(4) = {userdata.orientScore(userdata.idScore)-solution}; end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function snew = neigh_value_continuous(s,s_min,s_max,fraction1)

% select a RANDOM neighbouring real value of a SINGLE decision variable
% CEE 509, HW 5 by Bryan Tolson, Mar 5, 2003 AND ALSO CEE PROJECT

% variables:
% s is a current SINGLE decision variable VALUE
% s_min is the min of variable s
% s_max is the max of variable s
% snew is the neighboring VALUE of the decision variable
% fraction1 is the neighbourhood parameter (replaces V parameter~see notes)
%           It is defined as the ratio of the std deviation of the desired 
%           normal random number/s_range.  Eg:
%               std dev desired = fraction1 * s_range
%               for comparison:  variance (V) = (fraction1 * s_range)^2
% s_range is the range of the real variable (s_max-s_min)

s_range = s_max-s_min;

snew = s + randn(1,1) * fraction1 * s_range;

% NEED to deal with variable upper and lower bounds:
% Originally bounds in DDS were 100% reflective
% But some times DVs are right on the boundary and with 100% reflective
% boundaries it is hard to detect them. Therefore, we decided to make the
% boundaries reflective with 50% chance and absorptive with 50% chance.
% M. Asadzadeh and B. Tolson Dec 2008

P_Abs_or_Ref = rand();
if snew < s_min % works for any pos or neg s_min
    if P_Abs_or_Ref <= 0.5 %with 50% chance reflect
        snew = s_min + (s_min - snew); 
    else % with 50% chance absorb
        snew = s_min;
    end
        % if reflection goes past s_max then value should be s_min since without reflection
        % the approach goes way past lower bound.  This keeps X close to lower bound when X current
        % is close to lower bound:
    if snew > s_max
        snew = s_min;
    end
    
elseif snew > s_max % works for any pos or neg s_max
    if P_Abs_or_Ref <= 0.5%with 50% chance reflect
        snew = s_max - (snew - s_max); 
    else % with 50% chance absorb
        snew = s_max;
    end
        % if reflection goes past s_min then value should be s_max for same reasons as above
    if snew < s_min
        snew = s_max;
    end
end

% see 'neighbourhood bounding.xls' file for more analysis on above approach.

% Note the probability of sampling exactly at an endpoint based on St. Normal RV is:
%   = range because s current is fixed
%   = [Pr if at boundary of s Pr if at median of s]
%   = [Pmin(snew<min) + Pmin(snew>max), Pmax(snew<min) + Pmax(snew>max)]
%   = [    2*P(Z>2*s_range)           ,  2*P(Z>1.5*s_range)            ] %
%           Note: fraction1=stdev/s_range = 1/s_range, thus s_range = 1/fraction1
%   = [ 2*P(Z>2/fraction1), 2*P(Z>1.5/fraction1)]
% Table of resulting probabilities:
%   fraction1   Pmin        Pmax
%   <0.3        <0.0000     <0.0000 
%   0.4000      0.0000      0.0002
%   0.6000      0.0009      0.0124
%   0.8000      0.0124      0.0608
%   1.0000      0.0455      0.1336
%   1.2000      0.0956      0.2113

% above makes it difficult for neighbours to reach their absolute max and min 
% when fraction1 <0.3 but they can get close.  based on above, a
% fraction>1.2 does not seem to make very much sense since 1 of every 5 to 10
% samples of a parameter will result in visiting an exact endpoint.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function snew = neigh_value_discrete(s,s_min,s_max,fraction1)
% Created by B.Tolson and B.Yung, June 2006
% Modified by B. Tolson & M. Asadzadeh, Sept 2008
% Modification: 1- Boundary for reflection at (s_min-0.5) & (s_max+0.5)
%               2- Round the new value at the end of generation.
% select a RANDOM neighbouring integer value of a SINGLE decision variable
% discrete distribution is approximately normal 
% alternative to this appoach is reflecting triangular distribution (see Azadeh work)

% variables:
% s is a current SINGLE decision variable VALUE
% s_min is the min of variable s
% s_max is the max of variable s
% delta_s_min is the minimum perturbation size for each decision variable
    % equals [] if continuous DV (blank)
    % equals 1 if discrete integer valued DV
% snew is the neighboring VALUE of the decision variable
% fraction1 is the neighbourhood parameter (replaces V parameter~see notes)
%           It is defined as the ratio of the std deviation of the desired 
%           normal random number/s_range.  Eg:
%               std dev desired = fraction1 * s_range
%               for comparison:  variance (V) = (fraction1 * s_range)^2

% s_range is the range of the real variable (s_max-s_min)
s_range = s_max - s_min;
delta = randn(1,1) * fraction1 * s_range;
snew = s + delta;

P_Abs_or_Ref = rand();
if snew < s_min - 0.5 % works for any pos or neg s_min
    if P_Abs_or_Ref <= 0.5%with 50% chance reflect
        snew = (s_min-0.5) + ((s_min-0.5) - snew); 
    else %with 50% chance absorb
        snew = s_min;
    end
        % if reflection goes past (s_max+0.5) then value should be s_min since without reflection
        % the approach goes way past lower bound.  This keeps X close to lower bound when X current
        % is close to lower bound:
        if snew > s_max + 0.5
            snew = s_min;
        end        
elseif snew > s_max + 0.5 % works for any pos or neg s_max
    if P_Abs_or_Ref <= 0.5%with 50% chance reflect
        snew = (s_max+0.5) - (snew-(s_max+0.5));     
    else %with 50% chance absorb
        snew = s_max;
    end
        % if reflection goes past (s_min-0.5) then value should be s_max for same reasons as above
    if snew < s_min - 0.5
        snew = s_max;
    end
end

snew=round(snew); %New value must be integer
if snew==s %pick a number between s_max and s_min by a Uniform distribution
    sample = s_min-1+ceil((s_max-s_min)*rand(1,1)); % last term gives range = # options - 1.  First terms shift to allow min value
    if sample<s
        snew=sample;
    else  % must increment option number by one
        snew=sample+1;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function snew = neigh_value_mixed(s,s_min,s_max,fraction1,j,Discrete_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by B.Yung, April 2006
% Last modified by B.Yung, July 2006
% Purposes: run corresponding neighbour function for continuous and discrete 
%           decision variables  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% global Discrete_flag; %Discrete flag vector is entered in ext_function.inp

if Discrete_flag(j)==0  % continuous variable:
    snew=neigh_value_continuous(s,s_min,s_max,fraction1);
else                    % discrete integer variable
    snew=neigh_value_discrete(s,s_min,s_max,fraction1);    
end

end
